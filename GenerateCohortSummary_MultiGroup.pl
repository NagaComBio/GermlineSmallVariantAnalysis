#!/usr/bin/perl 

# Program to check case & control in HIPO_029
# Heterozygous variants

use strict;
use lib '/home/paramasi/perl5/lib/perl5';
use Config::Simple;
use Getopt::Long;

my %VARcmd;

# Get command line arguments
GetOptions ( "inputFilePrefix=s"       => \$VARcmd{'inputFilePrefix'},
             "inputFileSufix=s"       => \$VARcmd{'inputFileSufix'}, 
             "configFile=s"      => \$VARcmd{'configFile'},
             "outputFileSufixSummary=s"      => \$VARcmd{'outputFileSufix_Summary'},
             "outputFileSufixVP=s"      => \$VARcmd{'outputFileSufix_VP'},
             "baseDIR=s" => \$VARcmd{'baseDIR'},
             "project_ID=s" => \$VARcmd{'project_Id'}
           );

# Read the configuration file
my $cnf= new Config::Simple($VARcmd{'configFile'}) or die Config::Simple->error();

# Joining all the variants in HASH format
my %VARcnf = $cnf->vars();

my %VAR1 = (%VARcmd, %VARcnf);
my (%VAR) = CheckInputParameters(\%VAR1);
my $baseDIR = $VARcmd{'baseDIR'};

# Reading and writing the output files
open(IN, "<$baseDIR/$VAR{'inputFilePrefix'}$VAR{'inputFileSufix'}") || die "Cant open the $VAR{'inputFilePrefix'}$VAR{'inputFileSufix'}\n";
open(OUT, ">$baseDIR/$VAR{'inputFilePrefix'}$VAR{'outputFileSufix_Summary'}") || die "Cant create the $VAR{'inputFilePrefix'}$VAR{'outputFileSufix_Summary'}\n";

#### Basic info
my $basicInfoFile=$VAR{'basicInfo'};
open(BI, "$basicInfoFile") || die "cant open the basic info\n";
my %basicInfo;
while(<BI>)
{
  chomp;
  my @ss=split(/\t/, $_);
  #$ss[6]=~s/_X10//;
  $basicInfo{$ss[6]}{$ss[1]}{$ss[2]}++;
}

#######################################
my (%case, %caseFamGene, %caseGene, %caseALL);

my (%control, %controlFamGene, %controlGene, %controlALL);

my %gene;
my %annotation;
my %consensusVEP;

my %NuCaseWithAffectedGene;
my %NuControlWithAffectedGene;

my %geneFamilyAffected = 0;

my %caseUniqPID;
my %controlUniqPID;

my %caseUniqFamily;
my %controlUniqFamily;

my %caseUniqGender;
my %controlUniqGender;

my %caseUniqAge;
my %controlUniqAge;

my %VarTag_GM;
my %VarTagFamily_GM_CASE;
my %VarTagFamily_GM_CONTROL;
my %family_variants;

my @families;

my ($geneModel_Col, $CHR_Col, $POS_Col,$REF_Col, $ALT_Col, $PIDs_Col, $ANNOVAR_Col);
my ($FamilyID_Col, $GENE_Col, $EXONIC_Col, $var_Col, $Rare_Col, $Local_Col, $ESP_Col);
my ($metaSVM_col, $metaLR_col, $proveanPred_col, $CADD_col, $uniprotAcc_col);
my ($siftPred_col, $PPH2_HDIVPred_col, $PPH2_HVARPred_col, $LRT_Pred_col, $MT_Pred_col, $MA_Pred_col, $FATHMM_Pred_col);
my @dbNSFP_header;

my %allVariantsAllFamily;
my %controlSeen;

while(<IN>)
{
  chomp;
  my $line = $_;  
  if($_=~/^VAR_TYPE/)
  {
    my @splitHeader=split(/\t/, $line);
    for my $i(0 .. $#splitHeader)
    {
      $geneModel_Col = $i, if($splitHeader[$i]=~/^GeneticModel$/);
      $CHR_Col       = $i, if($splitHeader[$i]=~/^CHROM$/);
      $POS_Col       = $i, if($splitHeader[$i]=~/^POS$/);
      $REF_Col       = $i, if($splitHeader[$i]=~/^REF$/);
      $ALT_Col       = $i, if($splitHeader[$i]=~/^ALT$/);
      $PIDs_Col      = $i, if($splitHeader[$i]=~/^PatientPIDs$/);
      $FamilyID_Col  = $i, if($splitHeader[$i]=~/^Family_ID$/);
      $GENE_Col      = $i, if($splitHeader[$i]=~/^GENE$/);
      $var_Col       = $i, if($splitHeader[$i]=~/^VARIATION/);
      #$Rare_Col     = $i, if($splitHeader[$i]=~/^ExAC_Rareness/);
      $Local_Col     = $i, if($splitHeader[$i]=~/InAll_LocalControl-1-2-BlackList/);
      $ESP_Col       = $i, if($splitHeader[$i]=~/ESP6500/);
      
      # After this column all the information is either variant depend
      #  or gene depend, so it will be shared by all in the summary rows
      $EXONIC_Col = $i, if($splitHeader[$i]=~/^EXONIC_CLASSIFICATION$/);
      $ANNOVAR_Col = $i, if($splitHeader[$i]=~/^ANNOVARAnnotation$/);

      # dbNSFP columns
      $metaSVM_col     = $i, if($splitHeader[$i]=~/MetaSVM_pred/);
      $metaLR_col      = $i, if($splitHeader[$i]=~/MetaLR_pred/);
      $proveanPred_col = $i, if($splitHeader[$i]=~/PROVEAN_pred/);
      $CADD_col        = $i, if($splitHeader[$i]=~/CADD_PHRED_LOCAL/);
      $uniprotAcc_col  = $i, if($splitHeader[$i]=~/Uniprot_acc/);

      $siftPred_col      = $i, if($splitHeader[$i]=~/SIFT_pred/);
      $PPH2_HDIVPred_col = $i, if($splitHeader[$i]=~/Polyphen2_HDIV_pred/);
      $PPH2_HVARPred_col = $i, if($splitHeader[$i]=~/Polyphen2_HVAR_pred/);
      $LRT_Pred_col      = $i, if($splitHeader[$i]=~/LRT_pred/);
      $MT_Pred_col       = $i, if($splitHeader[$i]=~/MutationTaster_pred/);
      $MA_Pred_col       = $i, if($splitHeader[$i]=~/MutationAssessor_pred/);
      $FATHMM_Pred_col   = $i, if($splitHeader[$i]=~/FATHMM_pred/);

    }
    @dbNSFP_header = splice(@splitHeader, $EXONIC_Col);
  }
  else
  {
    print VP "$line\t";
    my @splitInputLine=split(/\t/, $line);

    $splitInputLine[0] =~s/var_//;
    my $phenotype = $splitInputLine[0];
    my $geneName = $splitInputLine[$GENE_Col];
    my $FamilyID = $splitInputLine[$FamilyID_Col];
    my $phenoFamilyID = $phenotype."_".$FamilyID;

    #Getting all the familyID in an array
    if(grep/^$phenoFamilyID$/, @families){}
    else{push(@families, $phenoFamilyID)}

    #Creating the VariantTag, Family specific variantTag, GeneticModel specific variantTag
    my $VariantTag = join("_", @splitInputLine[$CHR_Col..$ALT_Col]);
    my $FamVarTag = "$splitInputLine[$var_Col]"."_"."$phenoFamilyID"."#"."$VariantTag";   

    my $GM_FamVarTag = "$splitInputLine[$geneModel_Col]"."_"."$splitInputLine[$var_Col]"."_"."$FamVarTag";
    my $geneName_FamilyID = "$phenoFamilyID"."_"."$geneName";

    my $GM=$splitInputLine[$geneModel_Col];
    my $PID=$splitInputLine[$PIDs_Col];
   
    # Classification of variants based on Rareness in the population and local control
    my $Rareness=1;
    #if($splitInputLine[$Rare_Col]=~/RareAlleleExAC|^NA$/ && $splitInputLine[$Local_Col]=~/NoInBC/)
    #{$Rareness=1};
    
    #------------------------    
    $gene{$FamVarTag}=$geneName;
    
    # Adding all the annotation beyond exonic classification into hash 
    $annotation{$FamVarTag} = join("\t", @splitInputLine[$EXONIC_Col..$#splitInputLine]);
       
    # Family specific variant hash    
    $family_variants{$FamilyID}{$VariantTag}++;

    # Getting all the variant in list;
    $allVariantsAllFamily{$FamVarTag}++;
    
    push(@{$VarTag_GM{$VariantTag}}, $GM);
    #print "$FamVarTag\n"; 
    $case{$FamVarTag}++;
    $caseFamGene{$geneName_FamilyID}++;
    if($Rareness==1){$caseGene{$geneName}{$PID}{$GM}++;}
    $caseALL{$VariantTag}++;
    push(@{$caseUniqPID{$FamVarTag}}, $PID);
    push(@{$NuCaseWithAffectedGene{$geneName}}, $PID);  
    push(@{$VarTagFamily_GM_CASE{$FamVarTag}}, $GM);

    # Variant in Control phenotype
    if($phenotype =~/CONTROL/) {
      $controlSeen{$FamilyID}{$VariantTag}++;
    }
  }
}


# Printing the header for Cohort summary file
# New Columns


print OUT "Family_ID\tGENE\t";
print OUT "VARIATION\tCHROM_POS_REF_ALT\t";

foreach my $family (@families)
{ 
  print OUT "COUNT_", $family, "\t";
  print OUT "PID_", $family, "\t";
  print OUT "GT_", $family, "\t";
}

print OUT "InAllCase\tUniquenessOfVariant\t";

print OUT join("\t", @dbNSFP_header), "\n";
my %printed;

############# Based on Varaint ###################
foreach my $FamVarTag (keys %allVariantsAllFamily)
{
  my ($variation, $phenotype, $FamilyID, $VariantTag) = $FamVarTag=~/(\w+)_(\w+)_([Ff]amily.*)#(.*)/;  
  my ($CHR_POS) = $VariantTag=~/(\w\w?_\d+)_\w+_\w+/;
  my $newFamilyIDVariantTag = $FamilyID."_".$VariantTag;
  
  my $geneName = $gene{$FamVarTag};
  #$geneName=~s/\(.*\)//;

  my @uniqGM = uniq(@{$VarTag_GM{$VariantTag}});
  #my @uniqGM_Family = uniq(@{$VarTagFamily_GM{$FamVarTag}});
  my $uniqGM_tag = join("-", @uniqGM);
  my $uniqGM_Family_tag_CASE;

  # Family specific counts
  my $CASE_counts=0;
  if(defined $printed{$newFamilyIDVariantTag}){}
  else
  {
    print OUT "$FamilyID\t$geneName\t$variation\t$VariantTag\t";
  
    foreach my $otherFamily(@families)
    {
      my $new_FamVarTag=$variation."_".$otherFamily."#".$VariantTag;
 
      # Variant count in the family with a phenotype
      if(defined $case{$new_FamVarTag}){print OUT "$case{$new_FamVarTag}\t"}else{print OUT "0\t"}
    
      # Pids with the variants
      if(defined $caseUniqPID{$new_FamVarTag}){print OUT join(",", @{$caseUniqPID{$new_FamVarTag}}, "\t")} else{print OUT ".\t"}
      
      # Genotype of the variants
     if(defined $VarTagFamily_GM_CASE{$new_FamVarTag}){print OUT join(",", @{$VarTagFamily_GM_CASE{$new_FamVarTag}}, "\t")} else{print OUT ".\t"}
      
     ### Check if the variant is present in all the cases in the familiy
     #print "new\t$new_FamVarTag\n";
     if(defined $case{$new_FamVarTag})
     {
       my $caseFamilyTag="CASE_$FamilyID";
       #print "IN\t$FamVarTag\t$case{$new_FamVarTag}\t$basicInfo{$VARcmd{'project_Id'}}{$FamilyID}{'CASE'}\n";
       #print "$basicInfo{'CRC_X10'}{'Family_101'}{'CASE'}\t";
       #print "$basicInfo{'CRC_X10'}{'Family_102'}{'CASE'}\n";
       if($otherFamily=~/^$caseFamilyTag$/)
       {
         #print "INN\t$VARcmd{'project_Id'}\t$FamilyID\t";
         #print "$case{$new_FamVarTag}==$basicInfo{$VARcmd{'project_Id'}}{$FamilyID}{'CASE'}\n";
         if($case{$new_FamVarTag}==$basicInfo{$VARcmd{'project_Id'}}{$FamilyID}{'CASE'})
         {   
           $CASE_counts=1;
         }
         #elsif($case{$new_FamVarTag}==($basicInfo{$VARcmd{'project_Id'}}{$FamilyID}{'CONTROL'}))
         if(defined $controlSeen{$FamilyID}{$VariantTag})
         {
           $CASE_counts=2;
         }
       }
     }else{}
    }

    ## Variants in all cases or not  
    #print "$FamVarTag\t$CASE_counts\n";
    if($CASE_counts==1)
    {
      print OUT "AllCase-Yes\t";
    }
    elsif($CASE_counts==2){print OUT "AllControl-Yes-$controlSeen{$FamilyID}{$VariantTag}\t"}
    else
    {
      print OUT "AllCase-No\t";
    }
    # Number of exact variants from all the samples
    # if(defined $caseALL{$VariantTag}){print OUT "$caseALL{$VariantTag}\t"}else{print OUT "0\t"}
    # if(defined $controlALL{$VariantTag}){print OUT "$controlALL{$VariantTag}\t"}else{print OUT "0\t"}
  
    # Number of variants in the gene from all the sample
    # if(defined $caseGene{$geneName}){print OUT "$caseGene{$geneName}\t"}else{print OUT "0\t"}
    # if(defined $controlGene{$geneName}){print OUT "$controlGene{$geneName}\t"}else{print OUT "0\t"}

    my @uniqcaseGene = uniq(@{$NuCaseWithAffectedGene{$geneName}});
    my @uniqcontrolGene = uniq(@{$NuControlWithAffectedGene{$geneName}});

    # Tags the variant, whether it is present only in the family or not
    my @inWhichFamily=();
    map{if(defined $family_variants{$_}{$VariantTag}){push(@inWhichFamily, "$_")}} @families; 
    if(scalar(@inWhichFamily)>1){print OUT "NotUniqueToFamily\t"}
    elsif((scalar(@inWhichFamily))==1){if(grep(/$FamilyID/, @inWhichFamily)){print OUT "UniqueToFamily\t"}}
    else{print OUT "NA\t"}
    
    print OUT "$annotation{$FamVarTag}\n";
  }
  $printed{$newFamilyIDVariantTag}++;
}


sub uniq 
{
  return keys %{{ map { $_ => 1 } @_ }};
}

sub checkTolerance
{
  my $score=$_[0];
  my $itTag="NA";
  if($score <0){$itTag="Negative"}
  elsif($score >0){$itTag="Positive"}
  elsif($score == 0 ){$itTag="Zero"}
  elsif($score <-2){$itTag="VeryNegative"}
  elsif($score >2){$itTag="VeryPositive"}
  else{$itTag="NA"}

  return($itTag);
}

sub CheckInputParameters
{
  my %VARtoCheck = %{$_[0]};
  my %VARout;	
  while (my ($key, $value)=each %VARtoCheck)
  {
    $key =~s/default.//;
    $VARout{$key}=$value
  }
  my @mandatoryParameters;
  map {defined $VARout{$_} ? "NA" : die USAGE()} @mandatoryParameters;
  return(%VARout);
}
