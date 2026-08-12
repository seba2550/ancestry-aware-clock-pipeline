process STRATIFIED_MODELS {
    label 'process_low'
    label 'container_r'

    publishDir "${params.outdir}/aim3_models", mode: 'copy'

    input:
    path analysis_df

    output:
    path "stratified_logistic_models.csv", emit: stratified_results

    script:
    """
    stratified_models.R \\
        --analysis-rds ${analysis_df} \\
        --output-csv stratified_logistic_models.csv
    """
}
