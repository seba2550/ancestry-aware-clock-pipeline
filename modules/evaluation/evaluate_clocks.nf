process EVALUATE_CLOCKS {
    label 'process_medium'
    label 'container_r'

    publishDir "${params.outdir}/evaluation", mode: 'copy'

    input:
    path coefs_dir
    path magenta_beta
    path magenta_meta

    output:
    path "magenta_all_predicted_ages_and_acceleration.csv", emit: predictions
    path "composition_evaluation_summary.csv"              , emit: summary

    script:
    """
    evaluate_clocks.R \\
        --coefs-dir ${coefs_dir} \\
        --magenta-beta ${magenta_beta} \\
        --magenta-meta ${magenta_meta} \\
        --output-predictions magenta_all_predicted_ages_and_acceleration.csv \\
        --output-summary composition_evaluation_summary.csv
    """
}
