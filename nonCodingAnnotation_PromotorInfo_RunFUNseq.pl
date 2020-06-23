#!/usr/bin/perl
# preparing file for non-coding annotation 
# awk '{if($6==5){print $0}}' patient_CombinedGermlineAnnotation_Cohort_ExAC_SUMMARY.txt  | grep -v CommonAlleleExAC | grep NoInBlackList | grep NoInControl_1 > patient_CombinedGermlineAnnotation_Cohort_ExAC_SUMMARY_PC5-NoCommon-NoBlack-NoLC1.txt

use strict;
use Getopt::Long;

my %VARcmd;

# Get command line arguments
GetOptions ( "inFile=s"          => \$VARcmd{'inFile'},
             "pipelineDIR=s"     => \$VARcmd{'pipelineDIR'},
             "baseDir=s"         => \$VARcmd{'RESULTS_DIR'},
             "chr=s"             => \$VARcmd{'chr'},
             "var=s"             => \$VARcmd{'var'}
           );

my $inFile=$VARcmd{'inFile'};
my $pipelineDIR=$VARcmd{'pipelineDIR'};
my $baseDIR=$VARcmd{'RESULTS_DIR'};
my $chr = $VARcmd{'chr'};
my $variation = $VARcmd{'var'};

my $RUN_ReadInputData = 1;
my $RUN_createBED = 1;
my $RUN_funseq = 1;
my $RUN_encode = 0;
my $RUN_mergeVcf = 1 ;
my $RUN_addMeth = 0 ;
my $RUN_filterResults = 1;

my %localScore;
my %localScoreExonic;
my %localScoreEncodeMotif;
my %localScoreMETH;

my %dataInput;

###---------------------------------------------------------------------
# Prepare the files 
# Make sure the txt file has header

my $variationNO;
my $chrNO;
my $posNO;
my $refNO;
my $altNO;

if($RUN_ReadInputData == 1 )
{
  print "Reading input file\n";
  open(IN, "<$baseDIR/$inFile.$chr.$variation.CADD.GTEx.dbNSFP.txt") || die "cant open summary txt file $baseDIR/$inFile.$chr.$variation.CADD.GTEx.dbNSFP.txt\n";
  while(<IN>)
  {
    chomp;
    my @ss=split(/\t/, $_);

    if($_=~/Family_ID/)
    {
      for(my $i=0;$i<=$#ss; $i++)
      {
        if($ss[$i]=~/VARIATION/){$variationNO=$i;}
        if($ss[$i]=~/CHROM/){$chrNO=$i;}
        if($ss[$i]=~/POS/){$posNO=$i;}
        if($ss[$i]=~/REF/){$refNO=$i;}
        if($ss[$i]=~/ALT/){$altNO=$i;}
      }
    }
    else
    {
      my @var=($ss[$chrNO], $ss[$posNO], $ss[$refNO], $ss[$altNO]);
      my $key = join("_", $ss[$chrNO], $ss[$posNO], $ss[$refNO], $ss[$altNO]);
   
      if($var[0] eq $chr && $variation eq $ss[$variationNO])
      {      
        $dataInput{$key}++;
      }
    }  
  }
}
close IN;

if($RUN_createBED == 1)
{
  print "Creating BED file\n";
  open(BED, ">$baseDIR/$inFile.$chr.$variation.bed") || die "cant create $baseDIR/$inFile.$chr.$variation.bed\n";

  foreach my $key (keys %dataInput)
  {
    my @var=split(/_/, $key);
    
    my $col2;
    my $refLength=length($var[2]);
    my $ref=$var[1]-1;
  
    if($var[2]!~/,/ && $var[3]!~/,/)
    {
      my $altLength=length($var[3]);
      my $alt=$ref+$refLength;
      if($refLength == $altLength)
      {
        print BED "chr$var[0]\t$ref\t$alt\t$var[2]\t$var[3]\n";
      }
      else
      {
        print BED "chr$var[0]\t$ref\t$alt\t-\t-\n";
      }
    }
    elsif($var[3]=~/,/)
    {
      my @alts=split(/,/, $var[3]);
      my $all_SNVs=0;
      foreach my $alt(@alts)
      {
        if(length($alt) == 1){}else{$all_SNVs++}
      }
     
      if($all_SNVs==0)
      {
        my $ref=$var[1]-1;
        print BED "chr$var[0]\t$ref\t$var[1]\t-\t-\n";  
      }
      else
      {
        my $alt=$var[1]+$refLength;
        print BED "chr$var[0]\t$ref\t$alt\t-\t-\n";
      } 
    }
    else
    {
      print "WARNING:Multi_ref\t$key\n";
    }
  }  
  close BED;
}

if($RUN_funseq == 1)
{
  # Running FunSeq
  print "Running FunSeq\n";
  
  system("$pipelineDIR/funseq2-1.2/funseq2.sh -f $baseDIR/$inFile.$chr.$variation.bed -maf 1 -m 2 -inf bed -outf vcf -db -nc -o $baseDIR/funseq2_out.$chr.$variation");
  
  `rm $baseDIR/$inFile.$chr.$variation.bed`
}

#-----------------------------------------------------------------------
# Data/File Merging 
# program to merge vcf output from funseq

if($RUN_mergeVcf == 1)
{
  print "Merging funseq results\n";
  open(VCF2, "<$baseDIR/funseq2_out.$chr.$variation/Output.vcf") || print "cant open $baseDIR/funseq2_out.$chr.$variation/Output.vcf\n";
  open(VCF4, "<$baseDIR/funseq2_out.$chr.$variation/Output.indel.vcf") || print "cant open $baseDIR/funseq2_out.$chr.$variation/Output.indel.vcf\n";
  
  open(OUT, ">$baseDIR/$inFile.$chr.$variation.CADD.GTEx.dbNSFP.FunSeq.txt") || die "cant create non coding annotation $inFile.$chr.$variation.CADD.GTEx.dbNSFP.FunSeq.txt";
 
  print "Reading funseq file\n"; 
  chomp(my @vcf2=<VCF2>);
  close VCF2;
  print "Reading funseq indel file\n";
  chomp(my @vcf4=<VCF4>);
  close VCF4;

  
  ############################
  # Reading the FUNSEQ file
  my %vcf2Data ;
  foreach my $funseq(@vcf2)
  {
    chomp($funseq);
    if($funseq!~/^#CHR/)
    {
      my @funSeqArray = split(/\t/, $funseq);
      
      $funSeqArray[0]=~s/chr//;
      $funSeqArray[1]=$funSeqArray[1];
      $vcf2Data{$funSeqArray[0]}{$funSeqArray[1]}=$funSeqArray[7];
    }
  }
  ############################
  ## Reading funseq indels
  foreach my $funseq_indel(@vcf4)
  {
    chomp($funseq_indel);
    if($funseq_indel !~/^#/)
    {
      my @funSeqArray = split(/\t/, $funseq_indel);

      $funSeqArray[0]=~s/chr//;
      $funSeqArray[1]=$funSeqArray[1];
      $vcf2Data{$funSeqArray[0]}{$funSeqArray[1]}=$funSeqArray[7]; 
    }
  }  

  #############################
  # Reading the SUMMARY file
  
  my $vcfHeader = "";

  open(IN, "<$baseDIR/$inFile.$chr.$variation.CADD.GTEx.dbNSFP.txt") || die "cant open summary txt file $inFile.$chr.$variation.CADD.GTEx.dbNSFP.txt\n";
  while(<IN>)
  {
    chomp;
    my $line = $_;
    if($line=~/Family_ID/)
    {
      print OUT "$line\tCodingVariantsOrNot";
      print OUT "\tNetworkHUB\tGeneUnderNegativeSelection\tNonCodingENCODEAnnotation";
      print OUT "\tMotifBreaking\tMotifGain\tInSensitiveRegion\tInUltra-SensitiveRegion";
      print OUT "\tTargetGene\tKnownCancerGene\tTF_highlyOccupiedRegion";
      print OUT "\tRecurrentInDB";
      print OUT "\tFunSEQ_NonCodingScore_exp\tFunSEQ_NonCodingScore_10f";
      print OUT "\tGERP\tNCENC_Elements\n";
    }
    else
    {
      my @ss=split(/\t/, $line);
      #my ($chr, $pos) = $ss[5]=~/(\d\d?|X|Y)_(\d+)/;
      my $chr = $ss[$chrNO];
      my $pos = $ss[$posNO];
      
      my $countFunSeqAnnotation = 0;
      my $NCDS_score=0;
      my $NCENC_seen=0;
 
      if(defined $vcf2Data{$chr}{$pos})
      {
        $countFunSeqAnnotation = 1;	
        my @funseqAnnotation = split(/;/, $vcf2Data{$chr}{$pos});
        if(grep(/NCENC/, @funseqAnnotation)){$NCENC_seen=1;}
 
        print OUT "$line";          
        print OUT "\t", toPrint("CDS", \@funseqAnnotation);
        #print OUT "\t", toPrint("VA", \@funseqAnnotation);
        print OUT "\t", toPrint("HUB", \@funseqAnnotation);
        print OUT "\t", toPrint("GNEG", \@funseqAnnotation);
 
        my $ncenc = toPrint("NCENC", \@funseqAnnotation);
        my $count_ncenc=()=$ncenc=~/,/g;
        if($NCENC_seen == 1){$count_ncenc++};
        print OUT "\t$ncenc";
 
        print OUT "\t", toPrint("MOTIFBR", \@funseqAnnotation);
        print OUT "\t", toPrint("MOTIFG", \@funseqAnnotation);
        print OUT "\t", toPrint("SEN", \@funseqAnnotation);
        print OUT "\t", toPrint("USEN", \@funseqAnnotation);
        print OUT "\t", toPrint("GENE", \@funseqAnnotation);
        print OUT "\t", toPrint("CANG", \@funseqAnnotation);
        print OUT "\t", toPrint("HOT", \@funseqAnnotation);
        #print OUT "\t", toPrint("RECUR", \@funseqAnnotation);
        print OUT "\t", toPrint("DBRECUR", \@funseqAnnotation);
        #print OUT "\t", toPrint("CDSS", \@funseqAnnotation);
 
        my $NCDS_score1=toPrint("NCDS", \@funseqAnnotation);
        $NCDS_score=sprintf("%.10f", $NCDS_score1);
        print OUT "\t$NCDS_score1\t$NCDS_score";
 
        print OUT "\t", toPrint("GERP", \@funseqAnnotation);
        print OUT "\t$count_ncenc\n";
      }
      elsif($countFunSeqAnnotation == 0)
      {
        print OUT $line;
        for(1..16){print OUT "\t.";}
        print OUT "\n";    
      }
    }
  }
  close OUT;
}

sub toPrint
{
  my $ID=$_[0];
  my @funseqArray=@{$_[1]};
  my $returnVariable ="";

  map{/^$ID/ ?  $returnVariable = $_ : "" } @funseqArray;

  if($returnVariable=~/^$/)
  {
    $returnVariable="NA";
  }
  else
  {
    $returnVariable=~s/$ID=//;
  }
  return ($returnVariable);
}

