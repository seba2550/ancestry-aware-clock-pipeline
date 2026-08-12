#!/usr/bin/env python3
import os
import csv
import random

test_dir = os.path.dirname(os.path.abspath(__file__)) + "/test_data"
os.makedirs(test_dir, exist_ok=True)
assets_meqtl_dir = os.path.dirname(os.path.abspath(__file__)) + "/../assets/meqtl_lists"
os.makedirs(assets_meqtl_dir, exist_ok=True)
assets_prs_dir = os.path.dirname(os.path.abspath(__file__)) + "/../assets/prs_weights"
os.makedirs(assets_prs_dir, exist_ok=True)

random.seed(42)

# 1. CpG Probes
cpgs = [f"cg{i:08d}" for i in range(1, 501)]
with open(f"{test_dir}/synthetic_cpg_names.txt", "w") as f:
    f.write("\n".join(cpgs) + "\n")

# 2. EUR Metadata & Betas (25 samples)
eur_ids = [f"EUR_{i:03d}" for i in range(1, 26)]
with open(f"{test_dir}/synthetic_eur_meta.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["sample_id", "age", "ancestry"])
    for sid in eur_ids:
        writer.writerow([sid, round(random.uniform(20, 75), 1), "EUR"])

with open(f"{test_dir}/synthetic_eur_betas.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["CpG"] + eur_ids)
    for cpg in cpgs:
        writer.writerow([cpg] + [round(random.betavariate(2, 2), 4) for _ in eur_ids])

# 3. AFR Metadata & Betas (25 samples)
afr_ids = [f"AFR_{i:03d}" for i in range(1, 26)]
with open(f"{test_dir}/synthetic_afr_meta.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["sample_id", "age", "ancestry"])
    for sid in afr_ids:
        writer.writerow([sid, round(random.uniform(20, 75), 1), "AFR"])

with open(f"{test_dir}/synthetic_afr_betas.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["CpG"] + afr_ids)
    for cpg in cpgs:
        writer.writerow([cpg] + [round(random.betavariate(2, 2), 4) for _ in afr_ids])

# 4. MAGENTA Metadata & Betas (20 samples)
magenta_ids = [f"MAG_{i:03d}" for i in range(1, 21)]
cohorts = ["REAAADI", "PRADI", "PERUVIAN", "CuADI"]
with open(f"{test_dir}/synthetic_magenta_meta.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["IID", "Beta_ID", "AGE_OF_EXAM", "COHORT", "AD_status", "sex", "APOE_e4"])
    for sid in magenta_ids:
        writer.writerow([
            sid, sid, round(random.uniform(50, 85), 1),
            random.choice(cohorts), random.choice([0, 1]),
            random.choice(["M", "F"]), random.choice([0, 1, 2])
        ])

with open(f"{test_dir}/synthetic_magenta_betas.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["CpG"] + magenta_ids)
    for cpg in cpgs:
        writer.writerow([cpg] + [round(random.betavariate(2, 2), 4) for _ in magenta_ids])

# 5. Asset meQTL files
for fname, subset in [("cosmo_meqtls.csv", cpgs[:50]), ("genoa_meqtls.csv", cpgs[25:75]), ("epigen_meqtls.csv", cpgs[10:60])]:
    with open(f"{assets_meqtl_dir}/{fname}", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["cpg"])
        for c in subset:
            writer.writerow([c])

# 6. Asset PRS weight stub
with open(f"{assets_prs_dir}/PGS003958_hmPOS_GRCh37.txt", "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["rsID", "chr_name", "chr_position", "effect_allele", "other_allele", "effect_weight"])
    for i in range(1, 51):
        writer.writerow([f"rs{1000+i}", 1, 10000 + i*10, "A", "G", round(random.gauss(0, 0.1), 5)])

print("Synthetic test dataset generated successfully using standard library Python.")
