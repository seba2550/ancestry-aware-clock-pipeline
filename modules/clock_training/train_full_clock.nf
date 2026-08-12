process TRAIN_FULL_CLOCK {
    label 'process_medium'
    label 'container_python'

    publishDir "${params.outdir}/models", mode: 'copy'

    input:
    path binary_files

    output:
    path "full_clock_coef.csv", emit: coefs
    path "full_clock_summary.csv", emit: summary

    script:
    """
    unset PYTHONHOME PYTHONPATH
    train_full_clock.py \
        --input-prefix clock_combined_full \
        --alpha ${params.alpha} \
        --lambda-opt ${params.fixed_lambda} \
        --output-coefs full_clock_coef.csv \
        --output-summary full_clock_summary.csv
    """
}
