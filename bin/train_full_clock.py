#!/usr/bin/env python3
"""
==============================================================================
train_full_clock.py
Full Ancestry Combined Clock Training CLI (No train/test split)
==============================================================================
Trains a single ElasticNet model using ALL available pooled samples without
bootstrapping.
"""

import sys
import time
import argparse
import gc
import numpy as np
import pandas as pd
try:
    from glmnet import ElasticNet as GlmNetElasticNet
    HAS_GLMNET = True
except ImportError:
    from sklearn.linear_model import ElasticNet as SklearnElasticNet
    HAS_GLMNET = False
from pathlib import Path

def transform_age(ages, adult_age=20):
    return np.array([
        np.log((a + 1) / (adult_age + 1)) if a < adult_age
        else (a - adult_age) / (1 + adult_age)
        for a in ages
    ])

def anti_transform_age(exps, adult_age=20):
    return np.array([
        (1 + adult_age) * np.exp(e) - 1 if e < 0
        else (1 + adult_age) * e + adult_age
        for e in exps
    ])

def load_training_data(prefix):
    """Load training data from binary export files."""
    prefix_str = str(prefix)
    dims_file = f"{prefix_str}_X_dims.txt"
    bin_file = f"{prefix_str}_X.bin"
    cpg_file = f"{prefix_str}_cpg_names.txt"
    y_file = f"{prefix_str}_y.csv"

    for filepath in [dims_file, bin_file, cpg_file, y_file]:
        if not Path(filepath).is_file():
            raise FileNotFoundError(f"Input file not found: {filepath}")

    print(f"Loading data from prefix: {prefix_str}...")
    with open(dims_file, "r") as f:
        lines = f.read().strip().split("\n")
        n_samples, n_features = int(lines[0]), int(lines[1])

    print(f"Reading {n_samples} samples x {n_features} features...")
    X = np.fromfile(bin_file, dtype=np.float64).reshape(n_samples, n_features)

    with open(cpg_file, "r") as f:
        cpg_names = [line.strip() for line in f if line.strip()]

    y_df = pd.read_csv(y_file)
    y_raw = y_df["age_raw"].values

    return X, y_raw, cpg_names, y_df["sample_id"].values

def parse_args():
    parser = argparse.ArgumentParser(
        description="Train full combined ElasticNet clock model on all pooled samples."
    )
    parser.add_argument(
        "--input-prefix",
        type=str,
        required=True,
        help="Path prefix for binary data (e.g. /path/to/clock_combined_full)"
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.5,
        help="Float elastic net alpha (default: 0.5)"
    )
    parser.add_argument(
        "--lambda-opt",
        type=float,
        default=0.010474566337054866,
        help="Float regularization parameter (default: 0.010474566337054866)"
    )
    parser.add_argument(
        "--output-coefs",
        type=str,
        required=True,
        help="Path to save model coefficient CSV"
    )
    parser.add_argument(
        "--output-summary",
        type=str,
        required=True,
        help="Path to save summary CSV"
    )
    return parser.parse_args()

def main():
    try:
        args = parse_args()

        print("=" * 60)
        print("  Full Combined Clock Training (Fixed Lambda, single model)")
        print(f"  Lambda: {args.lambda_opt:.8f} | Alpha: {args.alpha}")
        print("=" * 60)

        start_total = time.time()

        # 1. Load Data
        X, y_raw, cpg_names, ids = load_training_data(args.input_prefix)

        print(f"  Memory allocated for X: {X.nbytes / 1e9:.2f} GB")
        print(f"  y_raw range: {y_raw.min():.1f} - {y_raw.max():.1f}")

        y_t = transform_age(y_raw)

        gc.collect()

        print("\n  Training ElasticNet with fixed lambda...")
        start_train = time.time()

        if HAS_GLMNET:
            model = GlmNetElasticNet(
                alpha=args.alpha,
                lambda_path=np.array([args.lambda_opt]),
                n_splits=0,
                standardize=True,
                fit_intercept=True,
                max_iter=100000,
                tol=1e-7,
                random_state=42
            )
            model.fit(X, y_t)
            coefs = model.coef_path_[:, 0]
            intercept = float(model.intercept_path_[0])
        else:
            model = SklearnElasticNet(
                alpha=args.lambda_opt,
                l1_ratio=args.alpha,
                fit_intercept=True,
                max_iter=100000,
                tol=1e-7,
                random_state=42
            )
            model.fit(X, y_t)
            coefs = model.coef_
            intercept = float(model.intercept_)
        
        train_time = time.time() - start_train
        n_nonzero = int(np.sum(coefs != 0))

        print(f"\n  Optimal Lambda: {args.lambda_opt:.6f}")
        print(f"  Non-zero CpGs:  {n_nonzero}")

        preds_t = X @ coefs + intercept
        preds_years = anti_transform_age(preds_t)
        mae = float(np.mean(np.abs(preds_years - y_raw)))
        r_matrix = np.corrcoef(preds_years, y_raw)
        r = float(r_matrix[0, 1]) if r_matrix.shape == (2, 2) else 0.0

        print(f"  In-sample MAE:  {mae:.2f} years")
        print(f"  In-sample r:    {r:.4f}")

        # 2. Save coefficients CSV
        mask = coefs != 0
        coef_df = pd.DataFrame({
            "cpg": ["(Intercept)"] + [cpg_names[i] for i in range(len(cpg_names)) if mask[i]],
            "coefficient": [intercept] + coefs[mask].tolist()
        })

        out_coefs_path = Path(args.output_coefs)
        out_coefs_path.parent.mkdir(parents=True, exist_ok=True)
        coef_df.to_csv(out_coefs_path, index=False)
        print(f"  Model coefficients saved to: {out_coefs_path}")

        # 3. Save summary CSV
        out_summary_path = Path(args.output_summary)
        out_summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_df = pd.DataFrame([{
            "model_name": "Combined_Full",
            "n_train": len(y_raw),
            "n_cpgs_nonzero": n_nonzero,
            "optimal_lambda": args.lambda_opt,
            "intercept": intercept,
            "insample_mae": mae,
            "insample_r": r,
            "training_time_sec": train_time
        }])
        summary_df.to_csv(out_summary_path, index=False)
        print(f"  Summary saved to: {out_summary_path}")

        print(f"\n{'='*60}")
        print(f"  Total Script Time: {time.time() - start_total:.1f}s")
        print(f"{'='*60}")

    except Exception as e:
        print(f"ERROR during full clock training: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
