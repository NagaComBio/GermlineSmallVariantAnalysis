#!/usr/bin/env perl
# Program to separate the results by family and version 
 

use strict;
use Getopt::Long;

my %VARcmd;

# Get command line arguments
GetOptions ( "inFile_VP=s" => \$VARcmd{'inFile_VP'},                                       
             "inFile_ComHetero=s" => \$VARcmd{'inFile_ComHetero'},
             "version=s" => \$VARcmd{'version'},
             "outputFileSufix_CH=s" => \$VARcmd{'outputFileSufix_CH'},
             "project=s" => \$VARcmd{'projectName'},
             "pipelineDIR=s" => \$VARcmd{'pipelineDir'},
             "baseDIR=s" => \$VARcmd{'baseDIR'}  
           );


my $CombinedFile  = "$VARcmd{'baseDIR'}/$VARcmd{'inFile_VP'}.txt";
my $CompHeterFile = "$VARcmd{'baseDIR'}/$VARcmd{'inFile_ComHetero'}.txt";
my $versionDIR = $VARcmd{'version'};
print "PEIOPE - $VARcmd{'pipelineDir'}\n";

my $outDir;
my @model;
my $comHeter;
my $filePRE;

my $genePanelFile;
my @genePanel;

open(ALL, "<$CombinedFile") || die "Cant open the combined file $CombinedFile\n";
open(CH, "<$CompHeterFile") || die "Cant open the comp hetero $CompHeterFile\n";
chomp(my @data=<ALL>);
chomp(my @data_CH=<CH>);

my $header ;
my %data_ALL;

#---------------------------------------
# to get patient ids
my %patientHash;
foreach my $line(@data)
{
  my @ss=split(/\t/, $line);
  if($ss[2]!~/PatientPIDs/){$patientHash{$ss[2]}=$ss[1];}  
}
my @patient=keys %patientHash;


#----------------------------------------
foreach my $pat(@patient)
{
  my $outDir = 	"$VARcmd{'baseDIR'}/$patientHash{$pat}/$VARcmd{'version'}";
  if(! -d $outDir){`mkdir -p $outDir`};
  	
  open(OUT_F, ">$outDir/$pat.filtered.txt") || die "Cant create the output file $outDir/$pat.filtered.txt";
  open(OUT_U, ">$outDir/$pat.unfiltered.txt") || die "Cant create the output file $outDir/$pat.unfiltered.txt";
  
  #------ Header --------------
  foreach my $line (@data)
  {
    if($line=~/^VAR_TYPE/)
    {
	  print OUT_F "PIDs_Gene\tVarTag_Parent\t$line\n"; 
	  print OUT_U "$line\n";
	}  
  }
    
  foreach my $dat(@data)
  {
    if(grep(/\t$pat\t/, $dat))
    { 		  	        
      if(grep(/\tHeterozygous\t/, $dat))
      {
        print OUT_U "$dat\n";
      }
      else
      {	
        print OUT_F "\t\t$dat\n";
        print OUT_U "$dat\n";
      }    
    }
  }

  foreach my $dat(@data_CH)
  {
    if(grep(/\t$pat\t/, $dat))
    {
	  print OUT_F "$dat\n";
	}
  }
  ### Creating excel output ##### 
  #`perl $VARcmd{'pipelineDir'}/csv2xlx.pl $outDir/$pat.filtered.txt $outDir/$pat.$VARcmd{'version'}.filtered.xls`;
  #`perl $VARcmd{'pipelineDir'}/csv2xlx.pl $outDir/$pat.unfiltered.txt $outDir/$pat.$VARcmd{'version'}.unfiltered.xls`;
}
