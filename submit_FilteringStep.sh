#!/bin/bash
set -x
#cd $PBS_SCRATCH_DIR/$PBS_JOBID

source $CONFIG_FILE
module load perl/5.20.2
perl $PIPELINE_DIR/FilteringRareGermlineVariants.pl \
            --configFile=$CONFIG_FILE \
            --inputFile=$input_vcf \
            --summaryFile=$summary_file\
            --patientPIDs=$PIDS \
            --variation=$VARIATION \
            --diseaseModel=$D_MODEL 
