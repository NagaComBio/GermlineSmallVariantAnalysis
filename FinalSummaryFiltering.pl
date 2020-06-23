#!/usr/bin/perl
use strict;
use Getopt::Long;

my %VARcmd;

# Get command line arguments
GetOptions ( "inFile=s"          => \$VARcmd{'inFile'},
             "pipelineDIR=s"     => \$VARcmd{'pipelineDIR'},
             "baseDir=s"         => \$VARcmd{'baseDIR'},
           );

my $baseDIR = $VARcmd{'baseDIR'};
my $inFile = $VARcmd{'inFile'};

## Cancer predisposition genes
open(CPG, "</desktop-home/paramasi/results/cancerPredispositionGenes_ALL.txt") || die "cant open the CPG file\n";
chomp(my @CPG=<CPG>);

## Filtering variants from high high classification
open(IN, "$baseDIR/$inFile.CADD.GTEx.dbNSFP.FunSeq.txt") || die "cant open $baseDIR/$inFile.CADD.GTEx.FunSeq.txt\n";
open(FILTER, ">$baseDIR/$inFile.CADD.GTEx.dbNSFP.FunSeq.FILTERED.txt") || die "cant create $baseDIR/$inFile.CADD.GTEx.FunSeq.FILTERED.txt\n";

print "Filtering variants\n";
my $variationNO;
my $gerpNO;
my $gerpNO_2;
my $funseqNCScoreNO;
my $ncencElementsNO;
my $exonicClassNO;
my $annovarNO;
my $caddNO;
my $exonicConsensusNO;
my $gtexGeneClassNO;
my $promoterGeneNO;
my ($geneModel_Col, $CHR_Col, $POS_Col,$REF_Col, $ALT_Col, $PIDs_Col, $ANNOVAR_Col);
my ($FamilyID_Col, $GENE_Col, $EXONIC_Col, $var_Col, $Rare_Col, $Local_Col, $ESP_Col);
my ($metaSVM_col, $metaLR_col, $proveanPred_col, $CADD_col, $uniprotAcc_col);
my ($siftPred_col, $PPH2_HDIVPred_col, $PPH2_HVARPred_col, $LRT_Pred_col, $MT_Pred_col, $MA_Pred_col, $FATHMM_Pred_col);
my (@gtexGeneExpNO, $gtexGeneExpTissuesNO, $exacpLIwordNO, $exacMisZNO);
my ($rareness_col, $clinvar_col);

my ($clinvar_col);

##### reading in the InFile...
while(<IN>)
{
  chomp;
  my $line = $_;
  my @ss=split("\t",$line);

  if($line=~/Family_ID/)
  {
    print FILTER "$line\t";
    print FILTER "Local_Consensus_Variant_Effect_Prediction_Meta\t";
    print FILTER "Local_Consensus_Variant_Effect_Prediction_HighLow\t";
    print FILTER "Local_Consensus_Variant_Effect_Prediction_DeleteriousCount\t";
    print FILTER "FILTERED_REASON\tTier\tLevel\n";

    for(my $i=0;$i<=$#ss;$i++)
    {
      $GENE_Col      = $i, if($ss[$i]=~/^GENE$/);
      #$Rare_Col     = $i, if($ss[$i]=~/^ExAC_Rareness/);
      $Local_Col     = $i, if($ss[$i]=~/InAll_LocalControl-1-2-BlackList/);
      $ESP_Col       = $i, if($ss[$i]=~/ESP6500/);

      # After this column all the information is either variant depend
      #  or gene depend, so it will be shared by all in the summary rows
      $EXONIC_Col = $i, if($ss[$i]=~/^EXONIC_CLASSIFICATION$/);
      $ANNOVAR_Col = $i, if($ss[$i]=~/^ANNOVARAnnotation$/);

      # dbNSFP columns
      $metaSVM_col     = $i, if($ss[$i]=~/MetaSVM_pred/);
      $metaLR_col      = $i, if($ss[$i]=~/MetaLR_pred/);
      $proveanPred_col = $i, if($ss[$i]=~/PROVEAN_pred/);
      $CADD_col        = $i, if($ss[$i]=~/CADD_PHRED_LOCAL/);
      $uniprotAcc_col  = $i, if($ss[$i]=~/Uniprot_acc/);

      $siftPred_col      = $i, if($ss[$i]=~/SIFT_pred/);
      $PPH2_HDIVPred_col = $i, if($ss[$i]=~/Polyphen2_HDIV_pred/);
      $PPH2_HVARPred_col = $i, if($ss[$i]=~/Polyphen2_HVAR_pred/);
      $LRT_Pred_col      = $i, if($ss[$i]=~/LRT_pred/);
      $MT_Pred_col       = $i, if($ss[$i]=~/MutationTaster_pred/);
      $MA_Pred_col       = $i, if($ss[$i]=~/MutationAssessor_pred/);
      $FATHMM_Pred_col   = $i, if($ss[$i]=~/FATHMM_pred/);

      $clinvar_col       = $i, if($ss[$i]=~/Clinvarness/);
      $rareness_col      = $i, if($ss[$i]=~/Rareness/);

      if($ss[$i]=~/VARIATION/){$variationNO=$i;}
      if($ss[$i]=~/^GERP$/){$gerpNO=$i;}
      if($ss[$i]=~/^GERP\+\+_RS$/){$gerpNO_2=$i;}      
      if($ss[$i]=~/FunSEQ_NonCodingScore_10f/){$funseqNCScoreNO=$i;}
      if($ss[$i]=~/NCENC_Elements/){$ncencElementsNO=$i;}
      if($ss[$i]=~/EXONIC_CLASSIFICATION/){$exonicClassNO=$i;}
      if($ss[$i]=~/ANNOVARAnnotation/){$annovarNO=$i;}
      if($ss[$i]=~/CADD_PHRED_LOCAL/){$caddNO=$i;}
      if($ss[$i]=~/PromoterGene/){$promoterGeneNO=$i;}
      if($ss[$i]=~/GTEx_geneExpressionIn/){push(@gtexGeneExpNO, $i);}      
      if($ss[$i]=~/GTEx_geneNumberOfTissuesWithHighExpression/){$gtexGeneExpTissuesNO=$i;}
      if($ss[$i]=~/ExAC_pLI_word/){$exacpLIwordNO=$i}
      if($ss[$i]=~/ExAC_mis_z/){$exacMisZNO=$i}
    }
  }
  else
  {     
    # Classification of vartiants based on functional and functional effect prediction
    my $exonicClassification = $ss[$EXONIC_Col]."_".$ss[$ANNOVAR_Col];

    if($exonicClassification=~/stop/)
    {
      $exonicClassification="VOS_Stop";
    }
    elsif($exonicClassification=~/^frameshift/)
    {
      $exonicClassification="VOS_FS"
    }
    elsif($exonicClassification=~/splicing/ && $exonicClassification!~/ncRNA/)
    {
      $exonicClassification="VOS_Splice";
    }
    elsif($ss[$metaSVM_col] =~ /D/ || $ss[$metaSVM_col] =~/D/)
    {
      $exonicClassification="VOS_Meta"
    }
    else
    {
      $exonicClassification="VUS";
    }

    ## Classification of Variants effect prediction to high or low
    my $VEP_HighLow = "LOW";
    my $VEP_count = 0;

    if($ss[$siftPred_col] =~/D/ || ($ss[$PPH2_HDIVPred_col] =~/D|P/ || $ss[$PPH2_HVARPred_col] =~/D|P/) || $ss[$LRT_Pred_col] =~/D/ || $ss[$MT_Pred_col] =~/A|D/ || $ss[$MA_Pred_col] =~/H|M/ || $ss[$FATHMM_Pred_col] =~/D/ || $ss[$proveanPred_col] =~/D/)
    {
      $VEP_HighLow="HIGH";
    }
    elsif($ss[$siftPred_col] =~/\./ && $ss[$PPH2_HDIVPred_col] =~/\./ && $ss[$PPH2_HVARPred_col] =~/\./ && $ss[$LRT_Pred_col] =~/\./ && $ss[$MT_Pred_col] =~/\./ && $ss[$MA_Pred_col] =~/\./ && $ss[$FATHMM_Pred_col] =~/\./ && $ss[$proveanPred_col] =~/\./)
    {
      $VEP_HighLow="No_Prediction";
    }

    if($exonicClassification=~/VOS_[FS|Stop|Splice]/)
    {
      $VEP_HighLow = "HIGH"
    }

    ## VEP counting
    if($ss[$siftPred_col] =~/D/){$VEP_count++;}
    if($ss[$PPH2_HDIVPred_col] =~/D|P/ || $ss[$PPH2_HVARPred_col] =~/D|P/){$VEP_count++;}
    if($ss[$LRT_Pred_col] =~/D/){$VEP_count++;}
    if($ss[$MT_Pred_col] =~/A|D/){$VEP_count++;}
    if($ss[$MA_Pred_col] =~/H|M/){$VEP_count++;}
    if($ss[$FATHMM_Pred_col] =~/D/){$VEP_count++;}
    if($ss[$proveanPred_col] =~/D/){$VEP_count++;}

    #------------------------    
 
    ### Filter reasons
    my $filter=0;
    my $reason="";

    ## CPGs
    my $geneName=$ss[$GENE_Col];
    $geneName=~s/\(.*\)//g;

    if($geneName=~/,/)
    {
      my @genes=split(/,/, $geneName);
      my $cpg=0;
      foreach my $gene(@genes)
      {
        print "$gene\n";
        $gene=~s/\(.*\)//g;
        print $gene, "\n";
        if(grep(/^$gene$/, @CPG)){$cpg++;}
      }
      if($cpg>0){$reason=$reason.";CPG"; $filter++;}
    }
    else
    {
      if(grep(/^$geneName$/, @CPG)){$reason=$reason.";CPG"; $filter++;}
    }

    ## funseq score
    if($ss[$funseqNCScoreNO] > 3)
    {
      $reason=$reason.";funseqNCScore_gt3";
      $filter++;
    }
    elsif($ss[$funseqNCScoreNO] > 1.5)
    {
      $reason=$reason.";funseqNCScore_gt1.5";
      $filter++;
    }
    elsif($ss[$funseqNCScoreNO] > 0)
    {
      $reason=$reason.";funseqNCScore_lt1.5";
      $filter++;
    }
    #### Gerp score
    if($ss[$gerpNO] > 4 || $ss[$gerpNO_2] > 4)
    {
      $reason=$reason.";GERP_gt4";
      $filter++;
    }
    elsif($ss[$gerpNO] > 2 || $ss[$gerpNO_2] > 2)
    {
      $reason=$reason.";GERP_gt2";
      $filter++;
    }
    elsif($ss[$gerpNO]=~/^\.$/ && $ss[$gerpNO_2] =~/^\.$/)
    {
      $reason=$reason.";NoGERP";
      $filter++;
    }
    else
    {
      $reason=$reason.";GERP_lt2";
      $filter++;
    }
    ### exonic and splice classification
    my $LOH=0;
    if($ss[$exonicClassNO]=~/stop[gain|loss]|^frameshift/ || $ss[$annovarNO]=~/splicing/i)
    {   
      $reason=$reason.";VOS";
      $LOH=1;
      $filter++;
    }
    ## CADD score
    if($ss[$caddNO] > 30)
    {
      $reason=$reason.";CADD_gt30";
      $filter++;
    }
    elsif($ss[$caddNO] > 20)
    {
      $reason=$reason.";CADD_gt20";
      $filter++;
    }
    elsif($ss[$caddNO] > 13)
    {
      $reason=$reason.";CADD_gt13";
      $filter++;
    }
    elsif($ss[$caddNO] > 5)
    {
      $reason=$reason.";CADD_lt13";
      $filter++;
    }
    elsif($ss[$caddNO] > 0)
    {
      $reason=$reason.";CADD_lt5";
      $filter++;
    }
    else
    {
      $reason=$reason.";NoCADD";
      $filter++;
    }
    ### dbNSFP prediction and nonsys
    if($VEP_HighLow=~/High/i && $ss[$exonicClassNO]=~/nonsynonymous/)
    {
      $reason=$reason.";VUS_ConsensusPredHigh";
      $filter++;
    }      
    ### non coding elements
    if($ss[$ncencElementsNO] > 5)
    {
      $reason=$reason.";NCelement_gt5";
      $filter++;
    }
    elsif($ss[$ncencElementsNO] > 1)
    {
      $reason=$reason.";NCelement_gt1";
      $filter++;
    }
    ### GTEx gene class
    my $tissueCount_Gtex=0;
    foreach my $gtexTissue(@gtexGeneExpNO)
    {
      if($ss[$gtexTissue]=~/HIGH/ && ($ss[$gtexGeneExpTissuesNO] < 10 && $ss[$gtexGeneExpTissuesNO] > 0))
      {
        $tissueCount_Gtex++;
      }
    }
    if($tissueCount_Gtex > 0)
    {
      $reason=$reason.";GTExHighExp_LessTissues";
      $filter++;
    }
   
    #### ExAC intolerance
    if($ss[$exacpLIwordNO]=~/Intolerant/ && $LOH==1) 
    {
      $reason=$reason.";ExAC.pLI-Intolerant";
      $filter++;
    }
    if($ss[$exacMisZNO] > 2 && $ss[$exonicClassNO]=~/nonsynonymous/)
    {
      $reason=$reason.";ExAC.misZ-Intolerant";
      $filter++;
    }
 
    ### In Promoter
    if($ss[$promoterGeneNO]!~/NotInPromoter/ && $ss[$ncencElementsNO] > 0)
    {
      $reason=$reason.";inPromWithNCelement:".$ss[$promoterGeneNO];
      $filter++;
    }
   
    if($ss[$rareness_col] =~ /Rare|NA/) {
       $reason=$reason.";RareVariant";
       $filter++;
    }
    else {
       $reason=$reason.";CommonVariant";
       $filter++;
    }
    ### clinvar 
    if($ss[$clinvar_col] !~/\./) {
      $reason=$reason.";ClinVARVariant";
      $filter++;
    }

    ### Tier and Levels
    my $tier = "Tier_0";
    my $level = "Level_0";

    if($LOH==1) {
      $tier = "Tier_1";
      $level = "Level_1";
      if($ss[$caddNO] >= 13) {
        $level = "Level_2";
        if($ss[$exacpLIwordNO]=~/Intolerant/) {
          $level = "Level_3";
        }
      }
    }
    elsif($ss[$exonicClassNO]=~/nonsynonymous/) {
      $tier = "Tier_2";
      $level = "Level_1";
      if($ss[$caddNO] >= 13) {
        $level = "Level_2";     
        if($ss[$exacMisZNO] > 2) {
          $level = "Level_3";
        }
      }
    }

    ### Printing
    if($filter > 0)
    {
      $reason=~s/;//;
      print FILTER "$line\t$exonicClassification\t$VEP_HighLow\t$VEP_count\t$reason\t$tier\t$level\n";
    }
    else
    {
      print FILTER "$line\t$exonicClassification\t$VEP_HighLow\t$VEP_count\tNotCandidate\t$tier\t$level\n";  
    }
  }
}
close FILTER;
