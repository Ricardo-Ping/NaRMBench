import os
import sys
import numpy as np
import torch
import joblib
import toml
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score, roc_curve, precision_recall_curve, auc
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
                        help='Modification type: choose from m1A, m7G, m6A, psU,m5C.',
                        required=False,default='m6A',choices=['m1A', 'm7G', 'm6A', 'psU','m5C'])
    return parser

def plot_roc_curve(y_true, y_score, save_dir, filename="test_roc_curve.pdf"):
    y_pred_avg = np.mean(y_score, axis=0)
    fpr, tpr, _ = roc_curve(y_true, y_pred_avg)
    roc_auc = auc(fpr, tpr)
    
    plt.figure()
    plt.plot(fpr, tpr, color='blue', lw=2, label=f'ROC AUC = {roc_auc:.2f}')
    plt.plot([0, 1], [0, 1], color='gray', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate', fontsize=14)
    plt.ylabel('True Positive Rate', fontsize=14)
    plt.legend(loc='lower right')
    
    plt.savefig(os.path.join(save_dir, filename))
    plt.close()

def plot_prauc_curve(y_true, y_score, save_dir, filename="test_prauc_curve.pdf"):
    y_pred_avg = np.mean(y_score, axis=0)
    precision, recall, _ = precision_recall_curve(y_true, y_pred_avg, pos_label=1)
    pr_auc = auc(recall, precision)
    
    plt.figure()
    plt.plot(recall, precision, color='blue', lw=2, label=f'PR AUC = {pr_auc:.2f}')
    #plt.plot([0, 1], [0, 1], color='gray', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall', fontsize=14)
    plt.ylabel('Precision', fontsize=14)
    plt.legend(loc='lower right')
    
    plt.savefig(os.path.join(save_dir, filename))
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
        os.makedirs(save_dir)

    train_info = dict()
    train_info["model_config"] = model_config
    train_info["train_config"] = train_config
    train_info["train_config"]["learning_rate"] = lr
    train_info["train_config"]["epochs"] = n_epoch
    train_info["train_config"]["save_per_epoch"] = save_per_epoch
    train_info["train_config"]["weight_decay"] = weight_decay
    train_info["train_config"]["number_of_validation_iterations"] = n_iterations
    train_info["train_config"]["seed"] = seed

    with open(os.path.join(save_dir, "train_info.toml"), 'w', encoding='utf-8') \
            as f:
        toml.dump(train_info, f)

    model = MILModel(model_config).to(device)
    train_dl, val_dl, test_dl = build_dataloader(train_config, args.n_processes,args.modification)

    train_ds = train_dl.dataset
    val_ds = val_dl.dataset
    test_ds = test_dl.dataset

    print("Train y counts:", np.bincount(train_ds.labels.astype(int)))
    print("Val   y counts:", np.bincount(val_ds.labels.astype(int)))
    print("Test  y counts:", np.bincount(test_ds.labels.astype(int)))


    optimizer = torch.optim.Adam(model.parameters(), lr=lr,
                                 weight_decay=weight_decay)

    criterion = build_loss_function(train_config['loss_function'])

    train_results, val_results = train(model, train_dl, val_dl, optimizer, n_epoch, device,
                                       criterion, save_dir=save_dir,
                                       save_per_epoch=save_per_epoch,
                                       n_iterations=n_iterations)

    joblib.dump(train_results, os.path.join(save_dir, "train_results.joblib"))
    joblib.dump(val_results, os.path.join(save_dir, "val_results.joblib"))

    selection_criteria = ['avg_loss', 'roc_auc', 'pr_auc']
    val_results = joblib.load(os.path.join(save_dir, "val_results.joblib"))
    for selection_criterion in selection_criteria:
        val_loss = [val_results[selection_criterion][i]
                    for i in range(0, len(val_results[selection_criterion]), save_per_epoch)]

        if selection_criterion in ('avg_loss', 'avg_loss_read'):
            best_model = (np.argmin(val_loss) + 1) * save_per_epoch
        else:
            best_model = (np.argmax(val_loss) + 1) * save_per_epoch

        state_dict = torch.load(os.path.join(save_dir, "model_states", str(best_model), "model_states.pt"))
        torch.save(state_dict, os.path.join(save_dir, "{}.pt".format(selection_criterion)))

        model.load_state_dict(state_dict)
        test_results = validate(model, test_dl, device, criterion, n_iterations)
        print("Criteria: {criteria} \t"
              "Compute time: {compute_time:.3f}".format(criteria=selection_criterion, compute_time=test_results["compute_time"]))
        print("Test Loss: {loss:.3f} \t"
              "Test ROC AUC: {roc_auc:.3f} \t "
              "Test PR AUC: {pr_auc:.3f}".format(loss=test_results["avg_loss"],
                                                 roc_auc=test_results["roc_auc"],
                                                 pr_auc=test_results["pr_auc"]))
        print("=====================================")
        joblib.dump(test_results, os.path.join(save_dir, "test_results_{}.joblib".format(selection_criterion)))
 
        all_labels = test_results['y_true']
        all_preds = test_results['y_pred']
         
                
        plot_roc_curve(all_labels, all_preds, save_dir, f"test_roc_curve_{selection_criterion}.pdf")
        plot_prauc_curve(all_labels, all_preds, save_dir, f"test_prauc_curve_{selection_criterion}.pdf")

if __name__ == "__main__":
    parser = argparser()
    # 允许直接解析命令行参数
    args = parser.parse_args()
    
    # 简单的运行前检查
    if not os.path.exists(args.train_config):
        print(f"错误: 找不到训练配置文件 {args.train_config}")
        sys.exit(1)
        
    main(args)