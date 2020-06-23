#!/bin/bash
# Simple bash script to call the cohort summary program
set -x

PREFIX_NAME=`echo ${cohortsCombinedFile%.txt}`
echo $PREFIX_NAME

source $CONFIG_FILE
module load perl/5.20.2
###############################################################################
# MultiGroup summary
perl $PIPELINE_DIR/GenerateCohortSummary_MultiGroup.pl \
 --inputFilePrefix $inFile \
 --baseDIR $RESULTS_DIR \
 --inputFileSufix .CADD.GTEx.dbNSFP.FunSeq.FILTERED.txt \
 --outputFileSufixSummary .CADD.GTEx.dbNSFP.FunSeq.FILTERED.SUMMARY.txt \
 --configFile $CONFIG_FILE \
 --project_ID $PROJECT_ID

###############################################################################
