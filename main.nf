#!/usr/bin/env nextflow

/*
========================================================================================
    Ancestry-Aware Methylation Clock Pipeline (Nextflow DSL2)
========================================================================================
    Pipeline Stages:
    1. Data Preparation (EUR+AFR dataset merge and binary matrix export)
    2. Clock Training (Scatter-gather composition gradient training + full clock)
    3. CpG Feature Characterization & Composite Figure 5
    4. Clock Evaluation on MAGENTA Cohort
    5. Polygenic Risk Score Integration & Aim 3 Multivariable / Quadrant AD Modeling
========================================================================================
*/

nextflow.enable.dsl = 2

// Import Modules
include { EXPORT_TRAINING_DATA } from './modules/data_prep/export_training_data'
include { TRAIN_ELASTIC_NET }    from './modules/clock_training/train_elastic_net'
include { TRAIN_FULL_CLOCK }     from './modules/clock_training/train_full_clock'
include { COLLECT_TRAINING_RESULTS } from './modules/clock_training/collect_results'
include { CHARACTERIZE_CPGS }    from './modules/characterization/characterize_cpgs'
include { GENERATE_FIG5 }        from './modules/characterization/generate_fig5'
include { EVALUATE_CLOCKS }      from './modules/evaluation/evaluate_clocks'
include { BUILD_ANALYSIS_DF }    from './modules/integration/build_analysis_df'
include { LOGISTIC_MODELS }      from './modules/integration/logistic_models'
include { STRATIFIED_MODELS }    from './modules/integration/stratified_models'
include { ELASTIC_NET_METALEARNER } from './modules/integration/elastic_net_metalearner'
include { QUADRANT_ANALYSIS }    from './modules/integration/quadrant_analysis'

workflow {
    log.info """
========================================================================================
  ANCESTRY-AWARE METHYLATION CLOCK PIPELINE v${workflow.manifest.version}
========================================================================================
  Run Profile         : ${workflow.profile}
  Output Directory    : ${params.outdir}
  Bootstrap Iterations: ${params.n_bootstrap}
  Elastic Net Alpha   : ${params.alpha}
  Fixed Lambda        : ${params.fixed_lambda}
========================================================================================
"""
    if (params.is_test_run) {
        log.info "Executing in TEST mode using synthetic test datasets..."
        
        // In test mode, create stubs or use synthetic data channels
        eur_beta  = Channel.fromPath("${params.test_data_dir}/synthetic_eur_betas.csv")
        eur_meta  = Channel.fromPath("${params.test_data_dir}/synthetic_eur_meta.csv")
        afr_beta  = Channel.fromPath("${params.test_data_dir}/synthetic_afr_betas.csv")
        afr_meta  = Channel.fromPath("${params.test_data_dir}/synthetic_afr_meta.csv")
        mag_beta  = Channel.fromPath("${params.test_data_dir}/synthetic_magenta_betas.csv")
        mag_meta  = Channel.fromPath("${params.test_data_dir}/synthetic_magenta_meta.csv")
    } else {
        // Production paths provided via params or GEO fetch
        eur_beta  = Channel.fromPath("${params.eur_beta}")
        eur_meta  = Channel.fromPath("${params.eur_meta}")
        afr_beta  = Channel.fromPath("${params.afr_beta}")
        afr_meta  = Channel.fromPath("${params.afr_meta}")
        mag_beta  = Channel.fromPath("${params.magenta_beta}")
        mag_meta  = Channel.fromPath("${params.magenta_meta}")
    }

    // 1. Data Preparation
    EXPORT_TRAINING_DATA(eur_beta, eur_meta, afr_beta, afr_meta)

    // 2. Scatter-Gather Clock Training
    // Define gradient ratios: tuple(ratio, eur_frac)
    ratios = Channel.of(
        tuple('100_0', 1.00),
        tuple('75_25', 0.75),
        tuple('50_50', 0.50),
        tuple('25_75', 0.25),
        tuple('0_100', 0.00)
    )

    iterations = Channel.of(0..(params.n_bootstrap - 1))

    // Build combination channel: tuple(ratio, eur_frac, iter_idx, seed)
    training_inputs = ratios
        .combine(iterations)
        .map { ratio, eur_frac, iter_idx ->
            def seed = 1000 + (eur_frac * 100).toInteger() * 100 + iter_idx
            return tuple(ratio, eur_frac, iter_idx, seed)
        }
        .combine(EXPORT_TRAINING_DATA.out.binary_data.map { files -> [files] })

    TRAIN_ELASTIC_NET(training_inputs)

    // Collect all trained composition model coefficients and summaries
    COLLECT_TRAINING_RESULTS(
        TRAIN_ELASTIC_NET.out.coefs.collect(),
        TRAIN_ELASTIC_NET.out.summaries.collect()
    )

    // Train Full Combined Clock
    TRAIN_FULL_CLOCK(EXPORT_TRAINING_DATA.out.binary_data)

    // 3. CpG Feature Characterization
    CHARACTERIZE_CPGS(COLLECT_TRAINING_RESULTS.out.combined_coefs, COLLECT_TRAINING_RESULTS.out.combined_coefs)
    GENERATE_FIG5(COLLECT_TRAINING_RESULTS.out.combined_coefs, COLLECT_TRAINING_RESULTS.out.combined_coefs)

    // 4. MAGENTA Evaluation
    EVALUATE_CLOCKS(
        COLLECT_TRAINING_RESULTS.out.combined_coefs.map { file -> file.parent },
        mag_beta,
        mag_meta
    )

    // 5. Aim 3 Integration & Multivariable AD Risk Modeling
    v2_merged    = Channel.fromPath("${params.test_data_dir}/synthetic_magenta_meta.csv")
    prs_csv      = Channel.fromPath("${params.prs_weights}")
    ancestry_txt = Channel.fromPath("${params.test_data_dir}/synthetic_magenta_meta.csv")

    BUILD_ANALYSIS_DF(
        v2_merged,
        EVALUATE_CLOCKS.out.predictions,
        prs_csv,
        ancestry_txt
    )

    LOGISTIC_MODELS(BUILD_ANALYSIS_DF.out.analysis_df)
    STRATIFIED_MODELS(BUILD_ANALYSIS_DF.out.analysis_df)
    ELASTIC_NET_METALEARNER(BUILD_ANALYSIS_DF.out.analysis_df)
    QUADRANT_ANALYSIS(BUILD_ANALYSIS_DF.out.analysis_df)
}
