#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
from pathlib import Path

import numpy as np
import torch
import joblib
import toml
import matplotlib.pyplot as plt

from sklearn.metrics import (
    roc_curve, precision_recall_curve, auc,
    confusion_matrix, precision_score, recall_score, f1_score,
    matthews_corrcoef, balanced_accuracy_score
)
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

# --- 第一步：先修正路径，再导入 m6anet 模块 ---
current_dir = os.path.dirname(os.path.abspath(__file__))
# 假设目录结构是 root/m6anet/scripts/train.py，我们要把 root 加入路径
root_dir = os.path.dirname(os.path.dirname(current_dir))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

# --- 第二步：现在可以安全地从 m6anet 导入了 ---
try:
    from m6anet.utils.builder import build_dataloader, build_loss_function
    from m6anet.utils.constants import DEFAULT_MODEL_CONFIG
    from m6anet.utils.training_utils import train, validate
    from m6anet.model.model import MILModel
except ImportError as e:
    print(f"导入失败！请检查 root_dir 是否正确: {root_dir}")
    print(f"错误详情: {e}")
    sys.exit(1)


def argparser():
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        add_help=False
    )
    parser.add_argument("--model_config",
                        help='path to model config file.',
                        default=DEFAULT_MODEL_CONFIG)
    parser.add_argument("--train_config",
                        help='path to training config file.',
                        required=True)
    parser.add_argument("--save_dir",
                        help='directory to output training results.',
                        required=True)
    parser.add_argument("--figure_dir",
                        help='directory (under save_dir) to save ROC/PR plots and npz files.',
                        default="figure",
                        type=str)

    parser.add_argument("--device",
                        help='device to perform training with.',
                        default='cpu', type=str)
    parser.add_argument("--lr",
                        help='training learning rate.',
                        default=4e-4, type=float)
    parser.add_argument("--seed",
                        help='random seed for training.',
                        default=25, type=int)
    parser.add_argument("--epochs",
                        help='number of training epochs',
                        default=50, type=int)
    parser.add_argument("--n_processes",
                        help='number of processes to use for training.',
                        default=25, type=int)
    parser.add_argument("--save_per_epoch",
                        help='number of epoch multiple to save training checkpoint.',
                        default=10, type=int)
    parser.add_argument("--weight_decay",
                        help='weight decay argument to regularize training.',
                        default=0, type=float)
    parser.add_argument("--num_iterations",
                        help='number of pass during evaluation step',
                        default=5, type=int)
    parser.add_argument('-m', '--modification',
                        help='Modification type: choose from m1A, m7G, m6A, psU, m5C.',
                        required=False, default='m6A',
                        choices=['m1A', 'm7G', 'm6A', 'psU', 'm5C'])
    return parser


# =========================
# 统一输出（对齐 TandemMod）
# =========================
def _to_pos_scores(y_score):
    """
    Convert y_score to 1D positive-class scores of shape (N,).
    Supports: (N,), (N,2), (T,N), (T,N,2)
    """
    y_score = np.asarray(y_score)

    if y_score.ndim == 1:  # (N,)
        return y_score

    if y_score.ndim == 2:
        # (N,2) or (T,N)
        if y_score.shape[1] == 2:      # (N,2)
            return y_score[:, 1]
        else:                           # (T,N)
            return y_score.mean(axis=0)

    if y_score.ndim == 3:
        # (T,N,2)
        y_mean = y_score.mean(axis=0)   # (N,2)
        if y_mean.shape[1] == 2:
            return y_mean[:, 1]
        raise ValueError(f"Unexpected y_score shape: {y_score.shape}")

    raise ValueError(f"Unsupported y_score shape: {y_score.shape}")


def export_curves_and_cm(
    y_true,
    y_score,
    out_dir,
    ckpt_base="m6A",
    title="m6A",
):
    """
    Save (same keys as your TandemMod script):
      {ckpt_base}.curves.npz : fpr, tpr, precision, recall
      {ckpt_base}.cm.npz     : tp, fp, fn, tn
    Plot (PDF):
      {title}_ROC_AUC.pdf
      {title}_PR.pdf
    """
    os.makedirs(out_dir, exist_ok=True)

    y_true = np.asarray(y_true).astype(int)
    pos_score = _to_pos_scores(y_score)  # (N,)

    # 用 0.5 阈值构造 y_pred（若你想用别的阈值可改这里）
    y_pred = (pos_score >= 0.5).astype(int)

    # Confusion matrix (force 2x2)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()

    # Curves
    fpr, tpr, _ = roc_curve(y_true, pos_score)
    roc_auc = auc(fpr, tpr)

    precision_curve, recall_curve, _ = precision_recall_curve(y_true, pos_score)
    pr_auc = auc(recall_curve, precision_curve)

    # Save npz (exact keys)
    npz_path = os.path.join(out_dir, f"{ckpt_base}.curves.npz")
    np.savez(npz_path, fpr=fpr, tpr=tpr, precision=precision_curve, recall=recall_curve)
    print(f"[saved] curves -> {npz_path}")

    cm_path = os.path.join(out_dir, f"{ckpt_base}.cm.npz")
    np.savez(cm_path, tp=tp, fp=fp, fn=fn, tn=tn)
    print(f"[saved] cm -> {cm_path}")

    # Print aligned metrics (optional)
    test_acc = (tp + tn) / max(1, (tp + tn + fp + fn))
    prec = precision_score(y_true, y_pred, zero_division=0)
    rec = recall_score(y_true, y_pred, zero_division=0)
    f1 = f1_score(y_true, y_pred, zero_division=0)
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    bal_acc = balanced_accuracy_score(y_true, y_pred)
    mcc = matthews_corrcoef(y_true, y_pred)

    print(
        "  Test :"
        f"  Accuracy {test_acc:.3f}  ROC-AUC {roc_auc:.3f}  PR-AUC {pr_auc:.3f}"
        f"  | Precision {prec:.3f}  Recall {rec:.3f}  F1 {f1:.3f}  MCC {mcc:.3f}"
        f"  | Specificity {specificity:.3f}  Balanced-Acc {bal_acc:.3f}"
    )

    # Plot ROC
    plt.figure()
    plt.plot(fpr, tpr, lw=2, label=f"AUC {roc_auc:.3f}")
    plt.plot([0, 1], [0, 1], linestyle="--", alpha=0.3)
    plt.xlabel("False positive")
    plt.ylabel("True positive")
    plt.title(title)
    plt.legend(loc="lower right")
    roc_pdf = os.path.join(out_dir, f"{title}_ROC_AUC.pdf")
    plt.savefig(roc_pdf, dpi=600, bbox_inches="tight")
    plt.close()

    # Plot PR
    plt.figure()
    plt.plot(recall_curve, precision_curve, lw=2, label=f"AUC {pr_auc:.3f}")
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title(title)
    plt.legend(loc="lower left")
    pr_pdf = os.path.join(out_dir, f"{title}_PR.pdf")
    plt.savefig(pr_pdf, dpi=600, bbox_inches="tight")
    plt.close()


def main(args):

    seed = args.seed
    np.random.seed(seed)
    torch.manual_seed(seed)

    device = args.device
    n_epoch = args.epochs
    lr = args.lr
    save_per_epoch = args.save_per_epoch
    save_dir = args.save_dir
    weight_decay = args.weight_decay
    n_iterations = args.num_iterations

    model_config = toml.load(args.model_config)
    train_config = toml.load(args.train_config)

    print("Saving training information to {}".format(save_dir))

    if not os.path.exists(save_dir):
        os.makedirs(save_dir, exist_ok=True)

    train_info = dict()
    train_info["model_config"] = model_config
    train_info["train_config"] = train_config
    train_info["train_config"]["learning_rate"] = lr
    train_info["train_config"]["epochs"] = n_epoch
    train_info["train_config"]["save_per_epoch"] = save_per_epoch
    train_info["train_config"]["weight_decay"] = weight_decay
    train_info["train_config"]["number_of_validation_iterations"] = n_iterations
    train_info["train_config"]["seed"] = seed

    with open(os.path.join(save_dir, "train_info.toml"), 'w', encoding='utf-8') as f:
        toml.dump(train_info, f)

    model = MILModel(model_config).to(device)
    train_dl, val_dl, test_dl = build_dataloader(train_config, args.n_processes, args.modification)

    train_ds = train_dl.dataset
    val_ds = val_dl.dataset
    test_ds = test_dl.dataset

    print("Train y counts:", np.bincount(train_ds.labels.astype(int)))
    print("Val   y counts:", np.bincount(val_ds.labels.astype(int)))
    print("Test  y counts:", np.bincount(test_ds.labels.astype(int)))

    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
    criterion = build_loss_function(train_config['loss_function'])

    train_results, val_results = train(
        model, train_dl, val_dl, optimizer, n_epoch, device,
        criterion,
        save_dir=save_dir,
        save_per_epoch=save_per_epoch,
        n_iterations=n_iterations
    )

    joblib.dump(train_results, os.path.join(save_dir, "train_results.joblib"))
    joblib.dump(val_results, os.path.join(save_dir, "val_results.joblib"))

    # ===== 选择最好 checkpoint 并测试 + 输出曲线 =====
    # selection_criteria = ['avg_loss', 'roc_auc', 'pr_auc']
    selection_criteria = ['roc_auc']
    val_results = joblib.load(os.path.join(save_dir, "val_results.joblib"))

    # 统一 figure 输出目录（对齐 TandemMod 的 figure_dir 概念）
    fig_dir = os.path.join(save_dir, args.figure_dir)
    os.makedirs(fig_dir, exist_ok=True)

    for selection_criterion in selection_criteria:
        val_loss = [
            val_results[selection_criterion][i]
            for i in range(0, len(val_results[selection_criterion]), save_per_epoch)
        ]

        if selection_criterion in ('avg_loss', 'avg_loss_read'):
            best_model = (np.argmin(val_loss) + 1) * save_per_epoch
        else:
            best_model = (np.argmax(val_loss) + 1) * save_per_epoch

        state_dict_path = os.path.join(save_dir, "model_states", str(best_model), "model_states.pt")
        state_dict = torch.load(state_dict_path, map_location=device)

        # 保存一个“最终选择”权重（你原来就这么做的）
        torch.save(state_dict, os.path.join(save_dir, "{}.pt".format(selection_criterion)))

        model.load_state_dict(state_dict)

        test_results = validate(model, test_dl, device, criterion, n_iterations)
        print(
            "Criteria: {criteria} \tCompute time: {compute_time:.3f}".format(
                criteria=selection_criterion,
                compute_time=test_results["compute_time"]
            )
        )
        print(
            "Test Loss: {loss:.3f} \tTest ROC AUC: {roc_auc:.3f} \tTest PR AUC: {pr_auc:.3f}".format(
                loss=test_results["avg_loss"],
                roc_auc=test_results["roc_auc"],
                pr_auc=test_results["pr_auc"]
            )
        )
        print("=====================================")

        joblib.dump(test_results, os.path.join(save_dir, f"test_results_{selection_criterion}.joblib"))

        # ====== 对齐 TandemMod 的曲线保存与 cm 保存 ======
        all_labels = test_results["y_true"]
        all_scores = test_results["y_pred"]  # 这里按“分数/概率/多次分数”处理

        # ckpt_base = f"{selection_criterion}"
        ckpt_base = "m6anet"
        export_curves_and_cm(
            y_true=all_labels,
            y_score=all_scores,
            out_dir=fig_dir,
            ckpt_base=ckpt_base,
            title=args.modification  # e.g., "m6A"
        )


if __name__ == "__main__":
    parser = argparser()
    args = parser.parse_args()

    # 简单的运行前检查
    if not os.path.exists(args.train_config):
        print(f"错误: 找不到训练配置文件 {args.train_config}")
        sys.exit(1)

    main(args)

"""
python /data/pingyc/projects/20250601_RNA/softwares/NaRMBench/2.retraining/2.1.retrain_m6Anet/m6anet/scripts/train2.py \
  --model_config /data/pingyc/projects/20250601_RNA/softwares/NaRMBench/2.retraining/2.1.retrain_m6Anet/m6anet/model/configs/model_configs/m6anet.toml \
  --train_config /data/pingyc/projects/20250601_RNA/softwares/NaRMBench/2.retraining/2.1.retrain_m6Anet/m6anet/model/configs/training_configs/m6anet_train_config.toml \
  --save_dir /data/pingyc/projects/20250601_RNA/data/HEK293T-S1/m6Anet/output \
  --device cuda \
  --lr 0.00001 \
  --seed 20 \
  --epochs 50 \
  --n_processes 25 \
  --save_per_epoch 10 \
  --num_iterations 1

"""
