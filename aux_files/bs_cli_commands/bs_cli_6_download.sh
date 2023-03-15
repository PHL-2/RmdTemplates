# download data
$HOME/bs.exe download project -i $1 \
    --exclude "*" \
    --include "*.zip" \
    --include "Undetermined_report_metrics.json" \
    --include "*detect.report.csv" \
    --include "*lineage_report.csv" \
    --include "*nextclade.tsv" \
    --include "*.consensus_hard_masked_sequence.fa" \
    --exclude "Undetermined.consensus_hard_masked_sequence.fa" \
    -o "$2" \
    --no-metadata \
    --overwrite