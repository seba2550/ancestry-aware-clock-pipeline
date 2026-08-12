process QUADRANT_ANALYSIS {
    label 'process_low'
    label 'container_r'

    publishDir "${params.outdir}/aim3_models", mode: 'copy'

    input:
    path analysis_df

    output:
    path "quadrant_concordance_stats.csv", emit: quadrant_stats

    script:
    """
    quadrant_analysis.R \\
        --analysis-rds ${analysis_df} \\
        --output-csv quadrant_concordance_stats.csv
    """
}
