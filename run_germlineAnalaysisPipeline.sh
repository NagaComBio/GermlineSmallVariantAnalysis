#!/bin/bash

## To run 
## One sample. 
#  sh run_germlineAnalaysisPipeline.sh -i PID
#
## For mulitple files
#  sh run_germlineAnalaysisPipeline.sh -i PID_1 PID_2 

  # Get the sample ID from command line
  
  while getopts 'vc:' OPTION
  do
      case $OPTION in
          v) vflag=1 # verbose
          ;;
	   c) iflag=1
	   CONFIG_FILE="$OPTARG"
	   ;;
          ?) printf "Usage: %s [-v] -i PATIENT_PID ...\n" $(basename $0) >&2
             exit 2
          ;;
      esac
  done
 
  shift $(($OPTIND - 1))
  CASE=$*
  
  # To keep the command line options simple, edit the config file path here
 # CONFIG_FILE="/home/paramasi/testing/configFiles/UniExome_config_germlineAnalaysisPipeline.sh"
  cd /
  ############### Testing / Checking
  # Checking for config file
  if [[ ! -f $CONFIG_FILE ]]
  then
      printf "Error: Config file %s not found. Please specify the absolute path to the config file.\n" $CONFIG_FILE >&2
      exit 2
  fi
  source $CONFIG_FILE

  # Checking for input PIDS 
  if [[  $CASE =~ ^$ ]]
  then
      printf "Error: The PID of the patient is not given.\n" >&2
      exit 2
  fi

  # Checking for pipeline DIR
  if [[ ! -d $PIPELINE_DIR ]]
  then
      printf "Analysis tools directory %s not found. Exiting...\n" $PIPELINE_DIR >&2
      exit 2
  fi
  
   
 
#####################################################################
# Starting the TRIO analysis pipeline
# Since SNVs and Indel are processed in parellel, the program has two
#  segments
####################################################################

CombinedFile_QSUB=()
CombinedFile_Names=()

if [[ $trioORsingle == 1 ]]
then
  FAMILY=patient
elif [[ $trioORsingle == 0 ]]   
then
  FAMILY=trio
elif [[ $trioORsingle == 2 ]] 
then
  FAMILY=cancer
fi

if [[ $RUN_FILTERING == 1 ]]
then
  for PIDS in ${CASE[@]}
  do
    echo "${PIDS}"
    Filter_QSUB=()
    Filter_Names=()
    
    ########################
    # Checking for folder to write result files
    eval "TrioAnalysisFolder=$TRIO_SUBDIR_PATH"
  
    if [ ! -d $TrioAnalysisFolder ]; then mkdir $TrioAnalysisFolder; fi
    cp $CONFIG_FILE $TrioAnalysisFolder
   
    ########################
    # Checkinf for the family results folder
    eval "ResultsPerFamily=$RESULTS_DIR" 
    if [ ! -d $ResultsPerFamily ]; then mkdir -p $ResultsPerFamily; fi 
     
    ######################
    # Filtering the variations ; i = D_model (disease model)
    for MODEL in ${RunMODEL[@]}
    do
      ####################
      # SNVs
      if [[ $RUN_SNV == 1 ]]
      then
        VARIATION=SNVs      
        if [[ $trioORsingle == 0 ]]
        then
          FAMILY=trio     
          eval "input_VCF=$P_SNV_VCF_FILE" # Input VCF file
          eval "summary_file=$SUMMARY_FILE" # Summary file
        fi
  
        if [[ $trioORsingle == 1 ]]
        then
          FAMILY=patient        
          eval "input_VCF=$P_SNV_VCF_FILE"
          eval "summary_file=$SUMMARY_FILE" # Summary file
        fi
  
        # For cancer samples. Platypus_SNV for Germline snvs and mpileup_snv for somatic_snvs
        if [[ $trioORsingle == 2 ]] 
        then
          FAMILY=cancer        
          if [[ $MODEL == "Somatic" ]] 
          then
            eval "input_VCF=$S_SNV_VCF_FILE" # somatic snvs
            eval "summary_file=$SUMMARY_FILE"            
            filter=`bsub -q medium -J SNV_FilterVCF_${MODEL}_${PIDS} -env "all,CONFIG_FILE=$CONFIG_FILE,PIDS=$PIDS,input_vcf=$input_VCF,D_MODEL=$MODEL,VARIATION=$VARIATION,summary_file=$summary_file" sh $PIPELINE_DIR/submit_FilteringStep.sh | cut -d'<' -f2 | cut -d'>' -f1`
            Filter_QSUB+=($filter)
            Filter_Names+=($summary_file)
            echo -e "\tFiltering for $MODEL $VARIATION : $filter"
          else
            eval "input_VCF=$P_SNV_VCF_FILE" # Platypus snv
            eval "summary_file=$SUMMARY_FILE"
          fi
        fi
        
        # Testing the input file
        if [[ ! -f $input_VCF ]]
        then
          printf "Error: Input SNVs VCF %s not found\n" $input_VCF>&2
          exit 2
        else     
          if [[ $MODEL != "Somatic" ]]     
          then            
            filter=`bsub -q medium -J SNV_FilterVCF_${MODEL}_${PIDS} -env "all,CONFIG_FILE=$CONFIG_FILE,PIDS=$PIDS,input_vcf=$input_VCF,D_MODEL=$MODEL,VARIATION=$VARIATION,summary_file=$summary_file" sh $PIPELINE_DIR/submit_FilteringStep.sh | cut -d'<' -f2 | cut -d'>' -f1`
            Filter_QSUB+=($filter)
            Filter_Names+=($summary_file)
            echo -e "\tFiltering for $MODEL $VARIATION : $filter"
          fi
        fi  
      fi
     
      #####################
      # INDELs       
      if [[ $RUN_INDEL == 1 ]]
      then
        VARIATION=Indel
        if [[ $trioORsingle == 0 ]]
        then
          FAMILY=trio
          eval "input_VCF=$INDEL_VCF_FILE" # Input VCF file
          eval "summary_file=$SUMMARY_FILE" # Summary file
        fi
  
        if [[ $trioORsingle == 1 ]]
        then
          FAMILY=patient        
          eval "input_VCF=$INDEL_VCF_FILE" # Input VCF file
          eval "summary_file=$SUMMARY_FILE" # Summary file                          
        fi
        
        if [[ $trioORsingle == 2 ]] 
        then
          FAMILY=cancer
          eval "summary_file=$SUMMARY_FILE"
          if [[ $MODEL == "Somatic" ]]        
          then
            eval "input_VCF=$S_INDEL_VCF_FILE" # Somatic Indels file
            eval "summary_file=$SUMMARY_FILE" # Summary file            
            filter=`bsub -q medium -J INDEL_FilterVCF_${MODEL}_${PIDS} -env "all,CONFIG_FILE=$CONFIG_FILE,PIDS=$PIDS,input_vcf=$input_VCF,D_MODEL=$MODEL,VARIATION=$VARIATION,summary_file=$summary_file" sh $PIPELINE_DIR/submit_FilteringStep.sh | cut -d'<' -f2 | cut -d'>' -f1`
            Filter_QSUB+=($filter)
            Filter_Names+=($summary_file)
            echo -e "\tFiltering for $MODEL $VARIATION : $filter"
          else
            eval "input_VCF=$INDEL_VCF_FILE" # For cancer germline indels
            eval "summary_file=$SUMMARY_FILE" # Summary file
          fi 
        fi
        # Testing the input file
        if [[ ! -f $input_VCF ]]
        then
          printf "Error: Input Indel VCF %s not found\n" $input_VCF>&2
          exit 2
        else  
          if [[ $MODEL != "Somatic" ]]
          then            
            filter=`bsub -q medium -J INDEL_FilterVCF_${MODEL}_${PIDS} -env "all,CONFIG_FILE=$CONFIG_FILE,PIDS=$PIDS,input_vcf=$input_VCF,D_MODEL=$MODEL,VARIATION=$VARIATION,summary_file=$summary_file" sh $PIPELINE_DIR/submit_FilteringStep.sh | cut -d'<' -f2 | cut -d'>' -f1`
            Filter_QSUB+=($filter)
            Filter_Names+=($summary_file)
            echo -e "\tFiltering for $MODEL $VARIATION : $filter"
          fi
        fi  
      fi
    done
    
  
    #############################
    ### Creating array for cohort summary
    eval "combinedFile=$COMBINED_OUTPUT_FILE"
    CombinedFile_Names+=($combinedFile)
  
  
    ###########################
    # Preparing combined file for a patient
    if [[ $RUN_join_GM_files == 1 ]]
    then      
      Filter_Jobs=$(printf "done(%s) &&" "${Filter_QSUB[@]}")      
      Filter_depend=`echo $Filter_Jobs | sed 's/ &&//'`      
      

      Filter_FileNames=$(printf "%s#" "${Filter_Names[@]}")

  
      joinFile=`bsub -q medium -J JoinFilterFile_${PIDS} -w "$Filter_depend" -env "all,Filter_FileNames=$Filter_FileNames,combinedFile=$combinedFile" sh $PIPELINE_DIR/joinFilteredFiles.sh | cut -d'<' -f2 | cut -d'>' -f1`
      CombinedFile_QSUB+=($joinFile)  
      echo -e "\tJoining disease models: $joinFile"
    fi  
  done 
fi

#################################
# Joining ALL pids
if [[ $RUN_JOINING_PIDS == 1 ]]
then  
  tmpFile=`shuf -i 10000-99999 -n 1`
  cm_tmpFile=${CLUSTER_EO}/${tmpFile}.txt
  (echo $CASE | tr " " "#") > $cm_tmpFile

  if [[ -z ${CombinedFile_QSUB[@]} ]]
  then
    CombinedFile_depend=""
  else    
    CombinedFile_Jobs=$(printf "done(%s) &&" "${CombinedFile_QSUB[@]}")
    CombinedFile_depend=`echo $CombinedFile_Jobs | sed 's/ &&$//'`    
  fi

  CombinedFile_FileNames=$(printf "%s#" "${CombinedFile_Names[@]}")
  eval "cohortsCombinedFile=$COHORT_COMBINED_FILE"  
  joiningPids=`bsub -q short -J ${PROJECT_ID}_joiningPids -w "$CombinedFile_depend" -env "all,CONFIG_FILE=$CONFIG_FILE,JOINED_PIDS=$cm_tmpFile,FAMILY=$FAMILY" sh $PIPELINE_DIR/submit_PidsConcatenation.sh | cut -d'<' -f 2 | cut -d'>' -f1`
  echo -e "Joining Pids: $joiningPids"
fi

eval "FINAL_SUMMARY=${FAMILY}_${TRIO_SUBDIR}_CombinedGermlineAnnotation_Cohort"
####################################
### Adding CADD, promotor, dbNSFP and FunSeq
if [[ $RUN_CADD_PROMOTER == 1 ]]
then

  if [ -z $joiningPids ]
  then
    joiningFile_depend=""
  else    
    joiningFile_depend="done($joiningPids)"
  fi

  echo -e "CADD and promoter jobs"
  CADD_PROMOTER_FunSeq_JOB=()
  
  for chr in `seq 1 22` X Y
  do
    for var in SNVs Indel
    do
      ## CADD 
      chr_var_file=${chr}'#'${var}'#'$FINAL_SUMMARY

      CADD_job_ID=`bsub -q long -J ${PROJECT_ID}_CADD_${chr}_${var} -w "$joiningFile_depend" -env "all,CONFIG_FILE=$CONFIG_FILE,chr_var_file=$chr_var_file" sh $PIPELINE_DIR/submit_CADDpromoter.run.sh | cut -d'<' -f 2 | cut -d'>' -f1`
      echo -e "\tCADD_Promoter & FunSeq job submitted for $var $chr: $CADD_job_ID"
      CADD_PROMOTER_FunSeq_JOB+=($CADD_job_ID)

    done
  done

  CADD_PROMOTER_FunSeq_JOBs=$(printf "done(%s) &&" "${CADD_PROMOTER_FunSeq_JOB[@]}")
  CADD_FunSeq_depend=`echo $CADD_PROMOTER_FunSeq_JOBs | sed 's/ &&$//'`


  ### Joining in  
  CADD_join_JOB=`bsub -q short -J ${PROJECT_ID}_CADDjoin_JOB -w "$CADD_FunSeq_depend" -env "all,PIPELINE_DIR=$PIPELINE_DIR,baseDIR=$RESULTS_DIR,inFile=$FINAL_SUMMARY" sh $PIPELINE_DIR/submit_CADDpromoter.JOIN.sh | cut -d'<' -f 2| cut -d'>' -f1`
  echo -e "CADD files join JOB: $CADD_join_JOB"
fi

##################################
# FILTERING 
FilteringJob=""
if [[ $RUN_FINAL_FILTERING == 1 ]]
then
  if [ -z $CADD_join_JOB ]
  then
    CombinedFile_depend=""
  else    
    CombinedFile_depend="done($CADD_join_JOB)"
  fi

  eval "cohortsCombinedFile=$COHORT_COMBINED_FILE"  
  FilteringJob=`bsub -q medium -J ${PROJECT_ID}_Filterin -w "$CombinedFile_depend" -env "all,PIPELINE_DIR=$PIPELINE_DIR,inFile=$FINAL_SUMMARY,baseDIR=$RESULTS_DIR" sh $PIPELINE_DIR/submit_FilteringJob.sh | cut -d'<' -f 2| cut -d'>' -f1`
  echo -e "Final Filtering : $FilteringJob"
fi

##################################
# Preparing summary on all Cohorts
CohortsSummary=""
if [[ $RUN_COHORT_ANALYSIS == 1 ]]
then

  if [ -z $FilteringJob ]
  then 
    FilteringJob_depend=""
  else    
    FilteringJob_depend="done($FilteringJob)"
  fi

  eval "cohortsCombinedFile=$COHORT_COMBINED_FILE"  
  CohortsSummary=`bsub -q medium -J ${PROJECT_ID}_CohortSummary -w "$FilteringJob_depend" -env "all,CONFIG_FILE=$CONFIG_FILE,inFile=$FINAL_SUMMARY" sh $PIPELINE_DIR/submit_CohortSummary.sh | cut -d'<' -f 2 | cut -d'>' -f1`
  echo -e "Cohort Summary: $CohortsSummary"
fi

####################################
# Generate Compound heterozygous and separete the results to filtered and unfilter
if [[ $RUN_CompHeter_Separate == 1 ]]
then

  if [ -z $CohortsSummary ]
  then 
    cohort_depend=""
  else    
    cohort_depend="done($CohortsSummary)"
  fi

  eval "cohortsCombinedFile=$COHORT_COMBINED_FILE"  
  CompoundHetero_Separate=`bsub -q medium -J ${PROJECT_ID}_CompHetero_Sep -w "$cohort_depend" -env "all,CONFIG_FILE=$CONFIG_FILE,inFile=$FINAL_SUMMARY" sh $PIPELINE_DIR/submit_compHetero_separate.sh | cut -d'<' -f 2 | cut -d'>' -f1`
  echo -e "CompoundHeterozygous and Separate: $CompoundHetero_Separate" 
fi
