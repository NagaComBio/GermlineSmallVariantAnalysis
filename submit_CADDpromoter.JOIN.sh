head $baseDIR/$inFile.1.SNVs.CADD.GTEx.dbNSFP.FunSeq.txt | grep Family_ID -m1 > $baseDIR/$inFile.CADD.GTEx.dbNSFP.FunSeq.txt

for chr in `seq 1 22` X Y
do
  for var in Indel SNVs
  do
    cat $baseDIR/$inFile.$chr.$var.CADD.GTEx.dbNSFP.FunSeq.txt | awk 'NR > 1' >> $baseDIR/$inFile.CADD.GTEx.dbNSFP.FunSeq.txt
    rm $baseDIR/$inFile.$chr.$var.CADD.GTEx.dbNSFP.FunSeq.txt
    rm $baseDIR/$inFile.$chr.$var.CADD.GTEx.dbNSFP.txt
    rm $baseDIR/$inFile.$chr.$var.CADD.GTEx.txt
    rm $baseDIR/$inFile.$chr.$var.CADD.txt
  done
done
