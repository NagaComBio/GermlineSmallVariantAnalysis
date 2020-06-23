#!/bin/bash
#
# Simple bash script to merge the files with differen genetic model

set -x
outputFileNames=($(echo $Filter_FileNames| tr "#" " "))

#if [[ -f $combinedFile ]]
#then
#  `rm $combinedFile`
#fi

# Overwring the $combinedFile

grep -h '^VAR_TYPE' ${outputFileNames[0]} > $combinedFile

for file in ${outputFileNames[@]}
  do    
#    grep -h '^VAR_C' $file >> $combinedFile        
    #cat $file | { grep '^VAR_C' || true; } >> "$combinedFile"
    (tail -n +2 $file || true ;) >> "$combinedFile"
    `rm $file`
  done  

