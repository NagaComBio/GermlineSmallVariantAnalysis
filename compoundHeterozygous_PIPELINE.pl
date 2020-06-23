#!/usr/bin/perl

# Program to get compound heterozygous 
# after adding functional impact annotation

use strict;
use Getopt::Long;

my %VARcmd;

# Get command line arguments
GetOptions ( "inFile=s" => \$VARcmd{'inFile'},             
             "outputFileSufix=s" => \$VARcmd{'outputFileSufix'},
             "trioORsingle=s" => \$VARcmd{'trioORsingle'},
             "baseDIR=s" => \$VARcmd{'baseDIR'}
           );

my $inFile  = "$VARcmd{'baseDIR'}/$VARcmd{'inFile'}.txt";
my $outFile = "$VARcmd{'baseDIR'}/$VARcmd{'inFile'}.$VARcmd{'outputFileSufix'}.txt";

open(IN, "<$inFile") || die "Cant open the heterozygous mutation file $inFile\n";
open(OUT, ">$outFile") || die "Cant create the Compound heterozygous mutations\n";

my %dataHash;
my %printHash;

my %dataHash_single;
my %geneCount_single;

my $headerPresent= 0;
my $header;
my ($chrNO, $posNO, $refNO, $altNO, $geneNO, $gmNO, $varNO, $rankNO);
my ($pidNO, $motherInfoNO, $fatherInfoNO);
my ($childGTNO, $motherGTNO, $fatherGTNO);

while(<IN>)
{
  chomp;
  my $line=$_;
  
  if($line=~/^VAR_TYPE/)
  {
    $headerPresent=1;
    $header=$line;

    my @ss=split(/\t/, $line);
    for(my $i=0;$i<=$#ss; $i++)
    {
      if($ss[$i]=~/CHROM/){$chrNO=$i;}
      if($ss[$i]=~/POS/){$posNO=$i;}
      if($ss[$i]=~/REF/){$refNO=$i;}
      if($ss[$i]=~/ALT/){$altNO=$i;}
      if($ss[$i]=~/GENE/){$geneNO=$i;}
      if($ss[$i]=~/GeneticModel/){$gmNO=$i;}
      if($ss[$i]=~/VARIATION/){$varNO=$i;}
      if($ss[$i]=~/Ranking/){$rankNO=$i;}
      if($ss[$i]=~/PatientPIDs/){$pidNO=$i;}      
      if($ss[$i]=~/ChildGT/){$childGTNO=$i;}
      if($ss[$i]=~/MotherGT/){$motherGTNO=$i;}
      if($ss[$i]=~/FatherGT/){$fatherGTNO=$i;}
    }
  }
  else
  {
    if($headerPresent==0)
    {
      die "Header missing in the Cohort file\n";  
    }
    my @split_input_line = split(/\t/, $line);        
   
    if($split_input_line[$gmNO]=~/Heterozygous/)
    {		
      my $geneName=$split_input_line[$geneNO];
      if($geneName=~/(\w+-?\w+)[;|,|\(]?/){$geneName=$1}
      
      my $PidGeneName = $split_input_line[$pidNO]."#".$geneName;
      my $VarTag=$split_input_line[$chrNO]."_".$split_input_line[$posNO]."_".$split_input_line[$refNO]."_".$split_input_line[$altNO];
  
      my $VarTagPid=join("_", $VarTag, $split_input_line[$pidNO]);
      
      $dataHash_single{$split_input_line[$pidNO]}{$geneName}{$VarTag}=$line;
      $geneCount_single{$split_input_line[$pidNO]}{$geneName}{$split_input_line[$motherGTNO]}++;      
            
      $printHash{$VarTagPid}=$line;
  
	    if($split_input_line[$motherGTNO]=~/HETER/ && $split_input_line[$fatherGTNO]=~/HOMO_R/)
	    {
	      my $value=$VarTag."_MOTHER";
	      push(@{$dataHash{$PidGeneName}},$value);
	    }
  
	    if($split_input_line[$motherGTNO]=~/HOMO_R/ && $split_input_line[$fatherGTNO]=~/HETER/)
	    {
	      my $value=$VarTag."_FATHER";
	      push(@{$dataHash{$PidGeneName}},$value);
	    }  
    } 
  }
}
close IN;

if($VARcmd{'trioORsingle'} == 0)
{
  while(my ($PidGeneName, $VarTags) = each %dataHash)
  {
    my $count ="";
    foreach my $VarTag(@$VarTags)
    {
      if($VarTag=~/MOTHER/) {$count="_MOTHER_".$count}
      if($VarTag=~/FATHER/) {$count="_FATHER_".$count}
    }
    if($count=~/MOTHER/ && $count=~/FATHER/)
    {
      foreach my $VarTag(@$VarTags)
      {
        my $VarTag_ID=$VarTag;      
      
        $VarTag =~ s/FATHER|MOTHER//g;
        my @subTag=split(/#/, $PidGeneName);
	      
        my $VarTagPid = join("", $VarTag, $subTag[0]);
        $printHash{$VarTagPid}=~s/Heterozygous/Compound Heterozygous/;
        print OUT "$PidGeneName\t$VarTag_ID\t$printHash{$VarTagPid}\n";
      }
    }
  }
}
elsif($VARcmd{'trioORsingle'} == 1)
{
  foreach my $pid (keys %geneCount_single)
  {
    foreach my $gene(keys %{$geneCount_single{$pid}})
    {
      if($geneCount_single{$pid}{$gene}{'HETER'}>=2)
      {
        foreach my $tag(keys %{$dataHash_single{$pid}{$gene}})
        {
          print OUT "$pid#$gene\t$tag\t$dataHash_single{$pid}{$gene}{$tag}\n";
          my $line = $dataHash_single{$pid}{$gene}{$tag};
          $line=~s/\tHeterozygous\t/\tCompoundHeterozygous\t/;
          print CM "$line\n";
        }
      }
      else
      {
        foreach my $tag(keys %{$dataHash_single{$pid}{$gene}})
        {
          print CM "$dataHash_single{$pid}{$gene}{$tag}\n";
        }
      }
    }
  }
}

close OUT;
