#!/bin/bash
## Config file for GermlineAanalysisPipeline
## To run, go to <tool_DIR>, sh run_germlineAnalaysisPipeline.sh PIDs

########################################################################
# ATTENTION:                                                           # 
#  The pipeline checks for the RARE, FUNCTIONAL variants in the VCF    # 
#  file. It will not distinguish between 'somatic' and 'germline'      #  
#  variants, so in the case of cancer samples, use the VCF file        # 
#  generated only from 'control/blood' sample.                         # 
########################################################################

########################################################################
# Variable that has to be parsed during runtime, will be given as      # 
# arguments to the perl script                                         #
########################################################################
#[UNSTABLE]
# Changes based on USER
PIPELINE_DIR=/home/$USER/projects/repository/germline-analysis-pipeline

## Root path for the data
DATA_ROOT_PATH=ROOT_PATH

# Changes with the corresponding project 
PROJECT_NAME=Dystonia_Pakistan
ANALYSIS_DIR=$DATA_ROOT_PATH/analysis/$PROJECT_NAME/sequencing/whole_genome_sequencing/results_per_pid
RESULTS_DIR=$DATA_ROOT_PATH/analysis/$PROJECT_NAME/sequencing/whole_genome_sequencing/results_per_family

# Analysis SubDir name
TRIO_SUBDIR=GA12_WES

# SubDir PWD
TRIO_SUBDIR_PATH=${ANALYSIS_DIR}/'${PIDS}'/${TRIO_SUBDIR}

### Cluster-related parameters
EMAIL=USER_EMAIL
CLUSTER_EO=$RESULTS_DIR/cluster_messages
ADD_CHILD_PBS_RESOURCES="walltime=00:30:00,nodes=1:ppn=1,mem=25g"
ADD_CADD_PBS_RESOURCES="walltime=05:00:00,nodes=1:ppn=1,mem=5g"
ADD_COHORT_PBS_RESOURCES="walltime=00:30:00,nodes=1:ppn=1,mem=5g"
ADD_JOINfiles_PBS_RESOURCES="walltime=00:10:00,nodes=1:ppn=1"
EO_MAIL_OPTS="-o $CLUSTER_EO -j oe -M $EMAIL -m a"

########################################################################
## Input files
# FAMILY will be 'trio' or 'patient' (based on trioORsingle)
# MODEL will be 'Denovo' 'Homozygous', 'Heterozygous' or 'Hemizygous' (based on RunMODEL)
# VARIATION will be 'SNVs' or 'Indel'
# 
# Program can take 'gzip'ped or 'gunzip'ped files

## TRIO input SNVs file
T_SNV_VCF_FILE=${ANALYSIS_DIR}/'${PIDS}'/addBAM/snvs_'${PIDS}'.vcf.gz_parentInfo.vcf

## Single patient input SNVs file
P_SNV_VCF_FILE=${ANALYSIS_DIR}/'${PIDS}'/mpileup/snvs_'${PIDS}'.vcf.gz

## TRIO or single Indel file
INDEL_VCF_FILE=${ANALYSIS_DIR}/'${PIDS}'/platypusIndel/indel'${PIDS}'.vcf.gz

## For Cancer Germline calls 
S_SNV_VCF_FILE=${ANALYSIS_DIR}/'${PIDS}'/${TRIO_SUBDIR}/platypus_snvs/snvs_'${PIDS}'.vcf.gz

########################################################################
## Output file
# if ==Summary, additional annotations will be added
# if ==VCF, only the variants passing the filters will be printed without any additional annotations
# The SUFFIX will be added in the main perl program
 
# Output Summary (==0) or VCF (==1)
toPrintSummaryORVCF=0
 
VCF_SUFFIX=_Filtered.vcf
SUMMARY_SUFFIX=_SummaryFiltered.txt

SUMMARY_FILE=${TRIO_SUBDIR_PATH}/'${FAMILY}'_'${MODEL}'_'${VARIATION}'_'${PIDS}'_SummaryFiltered.txt
COMBINED_OUTPUT_FILE=${TRIO_SUBDIR_PATH}/'${FAMILY}'_'${PIDS}'_CombinedGermlineAnnotation.txt
COHORT_COMBINED_FILE=${RESULTS_DIR}/'${FAMILY}'_'${TRIO_SUBDIR}'_CombinedGermlineAnnotation_Cohort.txt

########################################################################
## Run pipeline parts
# if ==0, Runs for trio (child and both parents)
# if ==1, Runs for single patients
# if ==2, Runs for single patients from cancer snv format
trioORsingle=1

# Variations SNVs and Indel
RUN_FILTERING=1

RUN_SNV=1
RUN_INDEL=1
RUN_join_GM_files=1

RUN_JOINING_PIDS=1
RUN_CADD_PROMOTER=1

RUN_COHORT_ANALYSIS=1
RUN_CompHeter_Separate=1

RUN_dbNSFP=1
RUN_nonCoding=1
RUN_MultiGroup=1

# Running with Confidence. Keep RUN_filterWithConCla as default for other
#  models, For cancer somatic, confidence and classification automatically 
#  applies

RUN_filterWithConCla=0
confidenceScore=7
CLASSIFICATION=somatic

# Available options for TRIO = Denovo, Homozygous, Heterozygous, Hemizygous
# Available options for Single = Homozygous, Heterozygous
# Available options for Cancer = Homozygous, Heterozygous, Somatic, LOH
# Use without ','
RunMODEL=(Homozygous Heterozygous)

## Disease Tissue
diseaseTissue='Brain;Nerve'

### Population Based Filters
runPopulationBasedFilter=0
subPopulation=''

########################################################################
# Stable arguments, will be used directly in Filtering_VCF.pl          #
########################################################################
#[ST]

# Base quality (QUAL)
QUAL=20

# Minimum read coverage (DP or DP4)
MinReads=10
ParentMinReads=10

# Minimum genotype quality score (GQ)
genotypeQuality=20

# Minor allele threshold for rare variant (AF in 1000Genome)
RareMAF=0.001
cutoff_LocalControl_2=0.1
cutoff_LocalControl_3=0.1
cutoff_LocalControl_4=0.1
cutoff_Phase3=0.01
cutoff_Phase1=0.01

# Variant allele frequency (Calculated from DP4 or DP5)
# Minimum VAF to be heterozygous
MinVAF_HT=0.10
# Maximum VAF to be heterozygous
MaxVAF_HT=0.90
# Maximum VAF to be Reference Homozygous
MaxVAF_HoRef=0.10
# Maximum VAF to be Alternative Homozygous
MaxVAF_HoAlt=0.90

# ANNOVAR annotations (patterns to be matched)
# SNP Exonic classification
patternSnpEc='nonsynonymous SNV|stop[gain|loss]'
# SNP ANNOVAR function # |\\bUTR
patternSnpAnn='splicing'
#patternSnpAnn='.'

# INDEL Exonic classification
patternIndelEc='frameshift [deletion|insertion|substitution]|stop[gain|loss]'
# INDEL ANNOVAR function |\\bUTR
patternIndelAnn='splicing'
#patternIndelAnn='.'

# Annotation with local Control and BlackList
# It requires considerable amount of time, so better to turn-off while testing.
# Prepare your own localControl, BlackList available.
RUN_LOCAL_CONTROL=0
RUN_BLACK_LIST=0

# Files used for additional annotations
# Basic info file, if it is TRIO, program gets the mother and father info from this file
# Use below format
#PID	FAMILY_ID	CASE_CONTROL	GENDER	MOTHER_ID	FATHER_ID	PROJECT
basicInfo=/home/$USER/results/basicInfo.txt

#Genetic assocation database
GADdb=$DATA_ROOT_PATH/analysis/UniExome/GermlineAnnotation/annotationInfo/genetic_association_DB/all.txt
#Gene description
geneDesc=$DATA_ROOT_PATH/analysis/UniExome/GermlineAnnotation/annotationInfo/proteinInfo/protein-coding_gene.txt
#OMIM disroders
omimGeneMap=$DATA_ROOT_PATH/analysis/UniExome/GermlineAnnotation/annotationInfo/OMIM/genemap2.txt
#GO annotations
GOannotation=GO_human/GO_from_MSigSB_forGenomeMusic.txt
#clinVar dbSNPs of pathogenic variants
clinVar=$DATA_ROOT_PATH/analysis/UniExome/GermlineAnnotation/annotationInfo/clinVar/clinvar_20160302.vcf
#DISEASE gene network
DisGeNET=$DATA_ROOT_PATH/analysis/UniExome/GermlineAnnotation/annotationInfo/DisGeNET/DisGeNET_CURATED_PARSED.txt

# intolerance scores 
RVIS_LocalControlFile=$DATA_ROOT_PATH/analysis/UniExome/GermlineAnnotation/annotationInfo/RVIS/RVIS_score_LocalControl.txt
RVIS_PlosGeneticsFile=$DATA_ROOT_PATH/analysis/UniExome/GermlineAnnotation/annotationInfo/RVIS/RVIS_score_PlosGenetics.txt
RVIS_ExAC=$DATA_ROOT_PATH/analysis/UniExome/GermlineAnnotation/annotationInfo/RVIS/RVIS_score_ExAC.r0.2.CalculatedFile.OnlyExonic.txt
GTEx_D=$DATA_ROOT_PATH/analysis/UniExome/GermlineAnnotation/annotationInfo/GTEx/GTEx_Analysis_V4_RNA-seq_RNA-SeQCv1.1.8_exon_reads.75Percentile.80PerSamples.txt.SUMMARY.bed

# Expression from bioGPS 
bioGPS_AnnFile=$DATA_ROOT_PATH/analysis/UniExome/GermlineAnnotation/annotationInfo/bioGPS/biogpsMAD_parsed.txt

# Subcellular localization and protein type
SCLinfoFile=$DATA_ROOT_PATH/analysis/UniExome/GermlineAnnotation/annotationInfo/proteinInfo/allGenes.Annotated.txt

# COSMIC mutations
COSMIC_CLEAN=$DATA_ROOT_PATH/analysis/UniExome/GermlineAnnotation/annotationInfo/COSMIC/CosmicMutantExport_v66_250713_MutationCoordinates_CLEAN.txt

## CADD Score
CADD_SNVs=$DATA_ROOT_PATH/analysis/UniExome/GermlineAnnotation/annotationInfo/CADD/whole_genome_SNVs.tsv.gz
CADD_Indel=$DATA_ROOT_PATH/analysis/UniExome/GermlineAnnotation/annotationInfo/CADD/InDels.tsv.gz

## Promoter
PromoterFile=$DATA_ROOT_PATH/analysis/UniExome/GermlineAnnotation/annotationInfo/promoter/gencode.v19.annotation.havana.promoter.gene_plain.bed.gz

## GTEX_tissueSpecific Expression file
gtexTissueSpecificExpFile=$DATA_ROOT_PATH/analysis/UniExome/GermlineAnnotation/annotationInfo/GTEx/GeneRPKM/GTEx_Analysis_V6_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.75Percentile.80PerSamples.txt
