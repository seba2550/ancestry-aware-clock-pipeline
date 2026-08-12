process EXPORT_TRAINING_DATA {
    label 'process_heavy_memory'
    label 'container_r'

    publishDir "${params.outdir}/data_prep", mode: 'copy'

    input:
    path eur_beta
    path eur_meta
    path afr_beta
    path afr_meta

    output:
    path "clock_combined_full*", emit: binary_data

    script:
    """
    export_training_data.R \\
        --eur-beta ${eur_beta} \\
        --eur-meta ${eur_meta} \\
        --afr-beta ${afr_beta} \\
        --afr-meta ${afr_meta} \\
        --output-prefix clock_combined_full
    """
}
