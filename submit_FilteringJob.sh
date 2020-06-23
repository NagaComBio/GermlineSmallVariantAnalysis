## Filtering
module load perl/5.20.2

perl $PIPELINE_DIR/FinalSummaryFiltering.pl \
  --baseDIR $baseDIR \
  --inFile $inFile \
  --pipelineDIR $PIPELINE_DIR
