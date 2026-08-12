process TRAIN_ELASTIC_NET {
    label 'process_medium'
    label 'container_python'

    input:
    tuple val(ratio), val(eur_frac), val(iter_idx), val(seed), path(binary_files)

    output:
    path "Composition_${ratio}_iter${iter_idx}_coef.csv", emit: coefs
    path "summary_${ratio}_iter${iter_idx}.csv", emit: summaries

    script:
    """
    train_elastic_net.py \\
        --input-prefix clock_combined_full \\
        --ratio ${ratio} \\
        --eur-frac ${eur_frac} \\
        --iteration ${iter_idx} \\
        --seed ${seed} \\
        --alpha ${params.alpha} \\
        --lambda-opt ${params.fixed_lambda} \\
        --output-coefs Composition_${ratio}_iter${iter_idx}_coef.csv \\
        --output-summary summary_${ratio}_iter${iter_idx}.csv
    """

    stub:
    """
    touch Composition_${ratio}_iter${iter_idx}_coef.csv
    touch summary_${ratio}_iter${iter_idx}.csv
    """
}
