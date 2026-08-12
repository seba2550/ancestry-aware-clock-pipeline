process LOGISTIC_MODELS {
    label 'process_low'
    label 'container_r'

    publishDir "${params.outdir}/aim3_models", mode: 'copy'

    input:
    path analysis_df

    output:
    path "primary_logistic_models.csv", emit: logistic_results

    script:
    """
    logistic_models.R \\
        --analysis-rds ${analysis_df} \\
        --output-csv primary_logistic_models.csv
    """
}
