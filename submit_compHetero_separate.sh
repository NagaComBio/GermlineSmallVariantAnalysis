#!/bin/bash
# Simple bash script to call compound heterozygous and separate the results
set -x
source $CONFIG_FILE
baseDIR=$RESULTS_DIR
module load perl/5.20.2
##############################################################################
# Generating compound heterozygous 

perl $PIPELINE_DIR/compoundHeterozygous_PIPELINE.pl \
 --inFile $inFile.CADD.GTEx.dbNSFP.FunSeq.FILTERED \
 --outputFileSufix ComHetero \
 --baseDIR $baseDIR \
 --trioORsingle $trioORsingle

#############################################################################
# Separating the files by family and induvidual pids

perl $PIPELINE_DIR/separateResults_PIPELINE.pl \
 --inFile_VP $inFile.CADD.GTEx.dbNSFP.FunSeq.FILTERED \
 --inFile_ComHetero $inFile.CADD.GTEx.dbNSFP.FunSeq.FILTERED.ComHetero \
 --baseDIR $baseDIR \
 --version $TRIO_SUBDIR \
 --project $PROJECT_ID \
 --pipelineDIR $PIPELINE_DIR 

#############################################################################
# Separate summary file by family 

perl $PIPELINE_DIR/separateSummaryBYFamily_PIPELINE.pl \
  --inFile $inFile.CADD.GTEx.dbNSFP.FunSeq.FILTERED.SUMMARY \
  --version $TRIO_SUBDIR \
  --outputDir $baseDIR \
  --pipelineDIR $PIPELINE_DIR
