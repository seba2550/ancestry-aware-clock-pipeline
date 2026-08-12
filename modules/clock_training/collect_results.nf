process COLLECT_TRAINING_RESULTS {
    label 'process_low'
    label 'container_python'

    publishDir "${params.outdir}/models", mode: 'copy'

    input:
    path coef_files
    path summary_files

    output:
    path "composition_v2_coefs.csv", emit: combined_coefs
    path "Composition_summary.csv", emit: combined_summary

    script:
    """
    unset PYTHONHOME PYTHONPATH
    python3 -c "
import pandas as pd
import glob

coef_files = sorted(glob.glob('Composition_*_iter*_coef.csv'))
dfs = []
for f in coef_files:
    ratio = f.split('_iter')[0].replace('Composition_', '')
    iteration = int(f.split('_iter')[1].replace('_coef.csv', ''))
    df = pd.read_csv(f)
    df['Ratio'] = ratio
    df['Iteration'] = iteration
    dfs.append(df)

if dfs:
    master_coefs = pd.concat(dfs, ignore_index=True)
    master_coefs.to_csv('composition_v2_coefs.csv', index=False)

sum_files = sorted(glob.glob('summary_*_iter*.csv'))
sum_dfs = [pd.read_csv(f) for f in sum_files]
if sum_dfs:
    master_sum = pd.concat(sum_dfs, ignore_index=True)
    master_sum.to_csv('Composition_summary.csv', index=False)
"
    """
}
