#!/usr/bin/perl 

use strict;
use Getopt::Long;

my %VARcmd;

# Get command line arguments
GetOptions ( "inFile=s"          => \$VARcmd{'inFile'},
             "pipelineDIR=s"     => \$VARcmd{'pipelineDIR'},
             "baseDir=s"         => \$VARcmd{'baseDir'},
             "chr=s"             => \$VARcmd{'chr'},
             "var=s"             => \$VARcmd{'var'},
             "GTEX_File=s"        => \$VARcmd{'gtexFile'},
             "diseaseTissue=s"     => \$VARcmd{'diseaseTissue'}
           );

my $chr = $VARcmd{'chr'};
my $variation = $VARcmd{'var'};
my $basePath="$VARcmd{'baseDir'}";
my $tissueString = $VARcmd{'diseaseTissue'};
my @tissueToSearch=split(/;/, $tissueString);

my $fileName="$VARcmd{'inFile'}.$chr.$variation.CADD";

########## - GTEx file
my $gtex="$VARcmd{'gtexFile'}";

open(GTEX, "<$gtex") || die "cant open the gtex file\n";

my %gtexExp;
my @TissueIndex;
while(<GTEX>)
{
  chomp;
  my @ss=split(/\t/, $_);
  if($ss[0]=~/^$/)
  {
    map{my $f=$_;push(@TissueIndex, grep $ss[$_] eq $f, 0 .. $#ss)} @tissueToSearch;
  }
  else
  {
    my @tissueExp;
    map{push(@tissueExp, $ss[$_])} @TissueIndex;
    for(my $i=0;$i<=$#tissueExp; $i++)
    {
      $gtexExp{$ss[0]}{$tissueToSearch[$i]}=$tissueExp[$i];
    }
    $gtexExp{$ss[0]}{'totalHigh'}=$ss[$#ss];
  }
}
close GTEX;

############### - ExAC files
my $exacFile = "ExAC/ExAC.r0.3.gene.mis_z.pLI.txt";
open(ExAC, "<$exacFile") || die "cant read in the $exacFile\n";
my %exacIntolerance;

while(<ExAC>)
{
  chomp;
  my @ss=split(/\t/, $_);
  my $decimal=sprintf("%.3f", $ss[2]);
  $exacIntolerance{$ss[0]}{'mis_z'}=$ss[1];
  $exacIntolerance{$ss[0]}{'pLI'}=$ss[2];

  if($decimal >= 0.9)
  {
    $exacIntolerance{$ss[0]}{'pLI_word'}="Intolerant";
  }
  else
  {
    $exacIntolerance{$ss[0]}{'pLI_word'}="Tolerant";
  }
}
close ExAC;

##### - reading in file
open(IN, "<$basePath/$fileName.txt") || die "cant open in file $basePath/$fileName.txt\n";
open(OUT, ">$basePath/$fileName.GTEx.txt") || die "cant create $basePath/$fileName.GTEx.txt\n";
my $geneIndex;

while(<IN>)
{
  chomp;
  my $line = $_;
  my @ss=split(/\t/, $line);
  
  if($line=~/Family_ID/)
  {
    ($geneIndex)  = grep $ss[$_] eq 'GENE', 0 .. $#ss;   
    print OUT "$line\t";
    foreach my $tissue(@tissueToSearch)
    {
      print OUT "GTEx_geneExpressionIn-$tissue\t";
    }
    print OUT "GTEx_geneNumberOfTissuesWithHighExpression\t";
    print OUT "ExAC_mis_z\tExAC_pLI\tExAC_pLI_word\n";
  }
  else
  {
    print OUT "$line\t";
    #### single gene
    my $geneName=$ss[$geneIndex];
    $geneName=~s/\(.*?\)//g;

    if($geneName!~/,/)
    {
      ### - GTEx expression
      if(defined $gtexExp{$geneName})
      {
        foreach my $tissue(@tissueToSearch)
        {
          print OUT "$gtexExp{$geneName}{$tissue}\t";
        }
        print OUT "$gtexExp{$geneName}{'totalHigh'}\t";

      }
      else
      {
        foreach my $tissue(@tissueToSearch)
        {
          print OUT "NoGTEx\t";
        }
        print OUT "NoGTEx\t";
      }
      
      ### ExAC intolerance
      if(defined $exacIntolerance{$geneName})
      {
        print OUT "$exacIntolerance{$geneName}{'mis_z'}\t";
        print OUT "$exacIntolerance{$geneName}{'pLI'}\t";
        print OUT "$exacIntolerance{$geneName}{'pLI_word'}\n";
      }
      else
      {
        print OUT ".\t.\t.\n";
      }
    }
    else
    {
      my @geneSplit = split(/,/, $geneName);
      ### GTEx
      foreach my $tissue(@tissueToSearch)
      {
        my @geneTissueExp;
        foreach my $gene(@geneSplit)
        { 
          if(defined $gtexExp{$gene})
          {
            push(@geneTissueExp, $gtexExp{$gene}{$tissue});
          }
          else
          {
            push(@geneTissueExp, "NoGTEx");
          }
        }
        print OUT join(",", @geneTissueExp), "\t";
      }
      ### Total tissues
      my $totalTissue="";
      foreach my $gene(@geneSplit)
      {
        if(defined $gtexExp{$gene})
        {
          $totalTissue="$totalTissue,$gtexExp{$gene}{'totalHigh'}";
        }
        else
        {
          $totalTissue="$totalTissue,NoGTEx";
        }
      }
      $totalTissue=~s/^,//;
      print OUT "$totalTissue\t";
 
      ### ExAC
      my @exacInto_pLI;
      my @exacInto_pLI_word;
      my @exacInto_mis_z;
      foreach my $gene(@geneSplit)
      {
        if($exacIntolerance{$geneName})
        {
          push(@exacInto_mis_z, $exacIntolerance{$geneName}{'mis_z'});
          push(@exacInto_pLI, $exacIntolerance{$geneName}{'pLI'});
          push(@exacInto_pLI_word, $exacIntolerance{$geneName}{'pLI_word'});
        }
        else
        {
          push(@exacInto_mis_z, ".");
          push(@exacInto_pLI, ".");
          push(@exacInto_pLI_word, ".");
        }
      }
      print OUT join(",", @exacInto_mis_z), "\t";
      print OUT join(",", @exacInto_pLI), "\t";
      print OUT join(",", @exacInto_pLI_word), "\n";
    }
  }
}
close IN;
close OUT;
