# Simple bash script to call the cohort summary program
set -x
source $CONFIG_FILE

###############################################################################
# Combine file
CombinedNames=($(cat $JOINED_PIDS | tr "#" " "))

#if [[ -f $cohortsSummaryFile ]]
#then
#  `rm $cohortsCombinedFile`
#fi

echo $cohortsCombinedFile
PIDS=${CombinedNames[0]}
eval "pidFile=$COMBINED_OUTPUT_FILE"
eval "cohortCombined=$COHORT_COMBINED_FILE"

#grep -h '^VAR_TYPE' ${CombinedFileNames[0]} | head -n 1 > $cohortsCombinedFile
grep -h '^VAR_TYPE' $pidFile | head -n 1 > $cohortCombined

for PIDS in ${CombinedNames[@]}
  do
    eval "pidFile=$COMBINED_OUTPUT_FILE"
    cat $pidFile | { grep '^var_' || true; } >> $cohortCombined
  done
