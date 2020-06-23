#!/bin/bash
set -x
source $CONFIG_FILE
module load htslib/1.4.1
module load bedtools/2.16.2
module load python/2.7.9
module load perl/5.20.2
module load cadd/1.3

chr_var_array=($(echo $chr_var_file| tr "#" " "))
chr=${chr_var_array[0]}
var=${chr_var_array[1]}
inFile=${chr_var_array[2]}

#eval inFile=$FINAL_SUMMARY

### CADD and promoter
perl $PIPELINE_DIR/CADD_Annotation.pl \
  --chr $chr \
  --var $var \
  --baseDIR $RESULTS_DIR \
  --inFile $inFile \
  --PIPELINE $PIPELINE_DIR \
  --promoterBED $PromoterFile

baseDIR=$RESULTS_DIR
#### Tissue specific annotations
perl $PIPELINE_DIR/AnnotateTissueSpecificExpression.pl \
  --chr $chr \
  --var $var \
  --baseDIR $baseDIR \
  --inFile $inFile \
  --pipelineDIR $PIPELINE_DIR \
  --GTEX_File $gtexTissueSpecificExpFile \
  --diseaseTissue $diseaseTissue

#### 
perl $PIPELINE_DIR/dbNSFP_Annotation.pl \
  --chr $chr \
  --var $var \
  --inputFilePrefix $inFile \
  --baseDIR $baseDIR \
  --pipelineDIR $PIPELINE_DIR \
  --outputFileSufix dbNSFP.txt \
  --populationBasedFilter $runPopulationBasedFilter \
  --ExACsubPopulation $subPopulation \
  --MAF $RareMAF

### FunSeq
perl $PIPELINE_DIR/nonCodingAnnotation_PromotorInfo_RunFUNseq.pl \
  --chr $chr \
  --var $var \
  --baseDIR $baseDIR \
  --inFile $inFile \
  --pipelineDIR $PIPELINE_DIR


### Filtering
#perl $PIPELINE_DIR/FinalSummaryFiltering.pl \
#  --chr $chr\
#  --var $var \
#  --baseDIR $baseDIR \
#  --inFile $inFile \
#  --pipelineDIR $PIPELINE_DIR
