#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pathlib
import warnings
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import numpy as np
import torch
import toml
from torch.utils.data import DataLoader

# --- 第一步：先修正路径，再导入 m6anet 模块 ---
current_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.dirname(os.path.dirname(current_dir))  # root/m6anet/scripts/inference2.py -> root
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

# --- 第二步：现在可以安全地从 m6anet 导入了 ---
try:
    from m6anet.model.model import MILModel
    from m6anet.utils.constants import (
        DEFAULT_MODEL_CONFIG,
        DEFAULT_MIN_READS,
        DEFAULT_READ_THRESHOLD,
        DEFAULT_NORM_PATH,
        PRETRAINED_CONFIGS,
        DEFAULT_PRETRAINED_MODEL,
        DEFAULT_PRETRAINED_MODELS,
    )
    from m6anet.utils.data_utils import NanopolishDS, NanopolishReplicateDS, inference_collate
    from m6anet.utils.inference_utils import run_inference
except ImportError as e:
    print(f"导入失败！请检查 root_dir 是否正确: {root_dir}")
    print(f"错误详情: {e}")
    sys.exit(1)


def argparser():
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        add_help=False
    )

    # Required
    parser.add_argument("--input_dir", nargs="*",
                        help="directories containing data.info and data.json.",
                        required=True)
    parser.add_argument("--out_dir",
                        help="directory to output inference results.",
                        required=True)

    # Optional
    parser.add_argument("--pretrained_model",
                        help=f"pre-trained model available at m6anet. Options include {DEFAULT_PRETRAINED_MODELS}.",
                        default=DEFAULT_PRETRAINED_MODEL, type=str)
    parser.add_argument("--model_config",
                        help="path to model config file.",
                        default=DEFAULT_MODEL_CONFIG)
    parser.add_argument("--model_state_dict",
                        help="path to model weights (.pt). If set, overrides pretrained_model.",
                        default=None)

    parser.add_argument("--norm_path",
                        help="path to normalization factors file",
                        default=DEFAULT_NORM_PATH)

    parser.add_argument("--batch_size",
                        help="batch size for inference.",
                        default=16, type=int)
    parser.add_argument("--save_per_batch",
                        help="saving inference results every save_per_batch multiples.",
                        default=2, type=int)
    parser.add_argument("--n_processes",
                        help="number of processes to run.",
                        default=25, type=int)
    parser.add_argument("--num_iterations",
                        help="number of sampling run.",
                        default=1000, type=int)
    parser.add_argument("--device",
                        help="device to perform inference with.",
                        default="cpu", type=str)
    parser.add_argument("--seed",
                        help="random seed for sampling.",
                        default=0, type=int)
    parser.add_argument("--read_proba_threshold",
                        help="default probability threshold for a read to be considered modified.",
                        default=DEFAULT_READ_THRESHOLD, type=float)

    parser.add_argument('-m', '--modification',
                        help="Modification type: choose from m1A, m7G, m6A, psU, m5C.",
                        required=False, default='m6A',
                        choices=['m1A', 'm7G', 'm6A', 'psU', 'm5C'])

    return parser


def main(args):
    # ---- model_state_dict / pretrained_model resolution ----
    if args.model_state_dict is not None:
        warnings.warn("--model_state_dict is specified, overwriting default pretrained model weights")
        # 你自己提供权重时：norm_path / read_proba_threshold 不会自动改
        # 请你自己确保这些与训练/预训练配置一致
    else:
        if args.pretrained_model not in DEFAULT_PRETRAINED_MODELS:
            raise ValueError(
                f"Invalid pretrained model {args.pretrained_model}, must be one of {DEFAULT_PRETRAINED_MODELS}"
            )
        # PRETRAINED_CONFIGS: name -> (state_dict_path, read_threshold, norm_path)
        args.model_state_dict = PRETRAINED_CONFIGS[args.pretrained_model][0]
        args.read_proba_threshold = PRETRAINED_CONFIGS[args.pretrained_model][1]
        args.norm_path = PRETRAINED_CONFIGS[args.pretrained_model][2]

    # ---- seeds ----
    torch.manual_seed(args.seed)
    if "cuda" in str(args.device):
        torch.cuda.manual_seed_all(args.seed)
    np.random.seed(args.seed)

    # ---- build model ----
    model_cfg = toml.load(args.model_config)
    model = MILModel(model_cfg).to(args.device)

    state = torch.load(args.model_state_dict, map_location=torch.device(args.device))
    model.load_state_dict(state)
    model.eval()

    # ---- output dir & headers ----
    pathlib.Path(args.out_dir).mkdir(parents=True, exist_ok=True)

    site_path = os.path.join(args.out_dir, "data.site_proba.csv")
    indiv_path = os.path.join(args.out_dir, "data.indiv_proba.csv")

    with open(site_path, 'w', encoding='utf-8') as f:
        f.write('transcript_id,transcript_position,n_reads,probability_modified,kmer,mod_ratio\n')
    with open(indiv_path, 'w', encoding='utf-8') as g:
        g.write('transcript_id,transcript_position,read_index,probability_modified\n')

    # ---- dataset ----
    if len(args.input_dir) == 1:
        ds = NanopolishDS(
            args.input_dir[0],
            modification=args.modification,
            min_reads=DEFAULT_MIN_READS,
            norm_path=args.norm_path,
            mode='Inference'
        )
    else:
        ds = NanopolishReplicateDS(
            args.input_dir,
            modification=args.modification,
            min_reads=DEFAULT_MIN_READS,
            norm_path=args.norm_path,
            mode='Inference'
        )

    dl = DataLoader(
        ds,
        num_workers=args.n_processes,
        collate_fn=inference_collate,
        batch_size=args.batch_size,
        shuffle=False
    )

    # ---- run inference ----
    run_inference(model, dl, args)
    print(f"[done] inference results saved to: {args.out_dir}")


if __name__ == "__main__":
    parser = argparser()
    args = parser.parse_args()
    main(args)
