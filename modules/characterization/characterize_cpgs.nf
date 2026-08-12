process CHARACTERIZE_CPGS {
    label 'process_medium'
    label 'container_r'

    publishDir "${params.outdir}/characterization", mode: 'copy'

    input:
    path comp_coefs
    path learning_coefs, stageAs: 'learning_coefs.csv'

    output:
    path "cpg_characterization_*", emit: stats
    path "core_genes_custom.csv" , emit: core_genes

    script:
    """
    characterize_cpgs.R \\
        --comp-coefs ${comp_coefs} \\
        --learning-coefs learning_coefs.csv \\
        --output-dir .

    for f in cpg_*.csv cpg_*.rds gene_*.csv gene_*.rds annotated_*.csv delta_beta_*.csv; do
        if [ -f "\$f" ]; then
            cp "\$f" "cpg_characterization_\$f"
        fi
    done

    if [ -f "core_genes.csv" ]; then
        cp core_genes.csv core_genes_custom.csv
    fi
    """

    stub:
    """
    touch cpg_characterization_stats.csv
    touch core_genes_custom.csv
    """
}
