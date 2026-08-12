process ELASTIC_NET_METALEARNER {
    label 'process_medium'
    label 'container_r'

    publishDir "${params.outdir}/aim3_models", mode: 'copy'

    input:
    path analysis_df

    output:
    path "metalearner_feature_weights.csv", emit: metalearner_weights

    script:
    """
    elastic_net_metalearner.R \\
        --analysis-rds ${analysis_df} \\
        --output-weights metalearner_feature_weights.csv
    """

    stub:
    """
    touch metalearner_feature_weights.csv
    """
}
