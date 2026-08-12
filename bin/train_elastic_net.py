#!/usr/bin/env python3
"""
==============================================================================
train_elastic_net.py
Single-Iteration Ancestry Composition Sensitivity Model Training CLI
==============================================================================
Trains a single ElasticNet clock model for a specified ancestry ratio and
bootstrap iteration using binary preprocessed methylation data.
"""

import sys
import time
import argparse
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
    """Load combined training data from binary export files."""
    prefix_str = str(prefix)
    dims_file = f"{prefix_str}_X_dims.txt"
    bin_file = f"{prefix_str}_X.bin"
    cpg_file = f"{prefix_str}_cpg_names.txt"
    y_file = f"{prefix_str}_y.csv"

    for filepath in [dims_file, bin_file, cpg_file, y_file]:
        if not Path(filepath).is_file():
            raise FileNotFoundError(f"Input file not found: {filepath}")

    print(f"Loading data from prefix: {prefix_str}")
    with open(dims_file, "r") as f:
        lines = f.read().strip().split("\n")
        n_samples, n_features = int(lines[0]), int(lines[1])

    X = np.fromfile(bin_file, dtype=np.float64).reshape(n_samples, n_features)

    with open(cpg_file, "r") as f:
        cpg_names = [line.strip() for line in f if line.strip()]

    y_df = pd.read_csv(y_file)
    y_raw = y_df["age_raw"].values

    return X, y_raw, cpg_names, y_df["sample_id"].values

def train_single(X, y_raw, cpg_names, alpha, fixed_lambda):
    """Train a single ElasticNet model (no CV)."""
    y_t = transform_age(y_raw)
    if HAS_GLMNET:
        model = GlmNetElasticNet(
            alpha=alpha,
            lambda_path=np.array([fixed_lambda]),
            n_splits=0,
            standardize=True,
            fit_intercept=True,
            max_iter=100000,
            tol=1e-7
        )
        model.fit(X, y_t)
        coefs = model.coef_path_[:, 0]
        intercept = float(model.intercept_path_[0])
    else:
        model = SklearnElasticNet(
            alpha=fixed_lambda,
            l1_ratio=alpha,
            fit_intercept=True,
            max_iter=100000,
            tol=1e-7
        )
        model.fit(X, y_t)
        coefs = model.coef_
        intercept = float(model.intercept_)

    preds_t = X @ coefs + intercept
    preds_years = anti_transform_age(preds_t)
    mae = float(np.mean(np.abs(preds_years - y_raw)))
    r_matrix = np.corrcoef(preds_years, y_raw)
    r = float(r_matrix[0, 1]) if r_matrix.shape == (2, 2) else 0.0

    mask = coefs != 0
    coef_df = pd.DataFrame({
        "cpg": ["(Intercept)"] + [cpg_names[i] for i in range(len(cpg_names)) if mask[i]],
        "coefficient": [intercept] + coefs[mask].tolist()
    })

    return coef_df, n_nonzero if 'n_nonzero' in locals() else int(np.sum(mask)), intercept, mae, r

def parse_args():
    parser = argparse.ArgumentParser(
        description="Train single-iteration ElasticNet composition clock model."
    )
    parser.add_argument(
        "--input-prefix",
        type=str,
        required=True,
        help="Path prefix for binary input data (e.g. /path/to/clock_combined_full)"
    )
    parser.add_argument(
        "--ratio",
        type=str,
        required=True,
        help="Training ratio label (e.g. 50_50, 100_0, 75_25, 25_75, 0_100)"
    )
    parser.add_argument(
        "--eur-frac",
        type=float,
        required=True,
        help="Float fraction of European samples (e.g. 0.50)"
    )
    parser.add_argument(
        "--n-total",
        type=int,
        default=1014,
        help="Int total training sample size (default: 1014)"
    )
    parser.add_argument(
        "--n-eur-available",
        type=int,
        default=507,
        help="Int available EUR samples (default: 507)"
    )
    parser.add_argument(
        "--n-afr-available",
        type=int,
        default=1380,
        help="Int available AFR samples (default: 1380)"
    )
    parser.add_argument(
        "--iteration",
        type=int,
        required=True,
        help="Int iteration index (e.g. 0 to 49)"
    )
    parser.add_argument(
        "--seed",
        type=int,
        required=True,
        help="Int random seed"
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
        help="Path to save output coefficient CSV"
    )
    parser.add_argument(
        "--output-summary",
        type=str,
        required=True,
        help="Path to save single-iteration summary CSV"
    )
    return parser.parse_args()

def main():
    try:
        args = parse_args()

        print("=" * 60)
        print(f"  Ancestry Composition Sensitivity (N={args.n_total}, Ratio={args.ratio}, Iteration={args.iteration})")
        print(f"  Lambda: {args.lambda_opt:.8f} | Alpha: {args.alpha} | Seed: {args.seed}")
        print("=" * 60)

        # 1. Load full sample pool
        X_pool, y_raw_pool, cpg_names, ids = load_training_data(args.input_prefix)
        y_df = pd.read_csv(f"{args.input_prefix}_y.csv")

        if "ancestry" in y_df.columns:
            idx_eur_pool = np.where(y_df["ancestry"] == "EUR")[0]
            idx_afr_pool = np.where(y_df["ancestry"] == "AFR")[0]
            if len(idx_eur_pool) == 0 and len(idx_afr_pool) == 0:
                mid = len(y_raw_pool) // 2
                idx_eur_pool = np.arange(0, mid)
                idx_afr_pool = np.arange(mid, len(y_raw_pool))
        else:
            mid = min(args.n_eur_available, len(y_raw_pool) // 2)
            idx_eur_pool = np.arange(0, mid)
            idx_afr_pool = np.arange(mid, len(y_raw_pool))

        n_total_target = min(args.n_total, len(y_raw_pool))
        n_eur_needed = int(round(n_total_target * args.eur_frac))
        n_afr_needed = n_total_target - n_eur_needed

        # 2. Bootstrap sampling
        rng = np.random.RandomState(args.seed)
        boot_idx_eur = rng.choice(idx_eur_pool, size=n_eur_needed, replace=True) if n_eur_needed > 0 and len(idx_eur_pool) > 0 else np.array([], dtype=int)
        boot_idx_afr = rng.choice(idx_afr_pool, size=n_afr_needed, replace=True) if n_afr_needed > 0 and len(idx_afr_pool) > 0 else np.array([], dtype=int)
        boot_idx = np.concatenate([boot_idx_eur, boot_idx_afr]) if len(boot_idx_eur) > 0 or len(boot_idx_afr) > 0 else np.arange(len(y_raw_pool))

        X_boot = X_pool[boot_idx]
        y_boot = y_raw_pool[boot_idx]

        # 3. Train ElasticNet model
        start_train = time.time()
        coef_df, n_nz, intercept, mae, r = train_single(
            X_boot, y_boot, cpg_names, args.alpha, args.lambda_opt
        )
        elapsed = time.time() - start_train

        # 4. Save coefficient CSV
        out_coefs_path = Path(args.output_coefs)
        out_coefs_path.parent.mkdir(parents=True, exist_ok=True)
        coef_df.to_csv(out_coefs_path, index=False)
        print(f"  Coefficients saved to: {out_coefs_path}")

        # 5. Save single-iteration summary CSV
        out_summary_path = Path(args.output_summary)
        out_summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_df = pd.DataFrame([{
            "condition": f"Ratio_{args.ratio}",
            "ratio_label": args.ratio,
            "eur_fraction": args.eur_frac,
            "iteration": args.iteration,
            "seed": args.seed,
            "n_train": args.n_total,
            "n_eur": n_eur_needed,
            "n_afr": n_afr_needed,
            "n_cpgs_nonzero": n_nz,
            "intercept": intercept,
            "insample_mae": mae,
            "insample_r": r,
            "training_time_sec": elapsed
        }])
        summary_df.to_csv(out_summary_path, index=False)
        print(f"  Summary saved to: {out_summary_path}")
        print(f"  Completed iteration {args.iteration} for ratio {args.ratio} (CpGs={n_nz}, MAE={mae:.2f}y, r={r:.4f}) in {elapsed:.1f}s")

    except Exception as e:
        print(f"ERROR during model training: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
