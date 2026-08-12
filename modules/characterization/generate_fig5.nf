process GENERATE_FIG5 {
    label 'process_medium'
    label 'container_r'

    publishDir "${params.outdir}/figures", mode: 'copy'

    input:
    path comp_coefs
    path learning_coefs, stageAs: 'learning_coefs.csv'

    output:
    path "Figure5.pdf", emit: fig5_pdf
    path "Figure5.png", emit: fig5_png

    script:
    """
    generate_fig5.R \\
        --comp-coefs ${comp_coefs} \\
        --learning-coefs learning_coefs.csv \\
        --output-dir .
    """

    stub:
    """
    touch Figure5.pdf
    touch Figure5.png
    """
}
