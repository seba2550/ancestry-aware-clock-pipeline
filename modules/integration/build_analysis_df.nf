process BUILD_ANALYSIS_DF {
    label 'process_low'
    label 'container_r'

    input:
    path v2_merged, stageAs: 'v2_merged_input.csv'
    path predictions, stageAs: 'predictions_input.csv'
    path prs_csv, stageAs: 'prs_csv_input.csv'
    path ancestry_txt, stageAs: 'ancestry_txt_input.txt'

    output:
    path "analysis_ready.rds", emit: analysis_df

    script:
    """
    build_analysis_df.R \
        --v2-merged v2_merged_input.csv \
        --custom-predictions predictions_input.csv \
        --prs-csv prs_csv_input.csv \
        --ancestry-txt ancestry_txt_input.txt \
        --output-rds analysis_ready.rds
    """
}
