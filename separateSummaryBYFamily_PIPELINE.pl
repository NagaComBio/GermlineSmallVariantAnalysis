#!/usr/bin/env perl
# Program to separate the results by family and version 
 

use strict;
use Getopt::Long;

my %VARcmd;

# Get command line arguments
GetOptions ( "inFile=s" => \$VARcmd{'inFile'},                                       
             "outputDir=s" => \$VARcmd{'outputDir'},
             "version=s" => \$VARcmd{'version'},                          
             "pipelineDIR=s" => \$VARcmd{'pipelineDir'}
           );


my $SummaryFile  = "$VARcmd{'outputDir'}/$VARcmd{'inFile'}.txt";
my $versionDIR = $VARcmd{'version'};

########################
open(IN, "<$SummaryFile") || die "cant open the summary file $SummaryFile\n";

chomp(my @data = <IN>);

##########################
my %FamilyHash;
foreach my $line(@data)
{
  my @ss=split(/\t/, $line);
  if($ss[0]!~/Family_ID/){$FamilyHash{$ss[0]}++;}  
}
my @Families=keys %FamilyHash;

################## 
if(scalar(@Families) == 0)
{
  die "Summary file without family Ids";
}
#elsif(scalar(@Families) ==1 )
#{
#  my @outFile = split(/\//, $VARcmd{'inFile'});    
#  my $familySummaryFile="$VARcmd{'outputDir'}/$Families[0]/$VARcmd{'version'}/$outFile[$#outFile].$Families[0]";
#  `cp $SummaryFile $familySummaryFile.txt`;
#  #`perl $VARcmd{'pipelineDir'}/csv2xlx.pl $familySummaryFile.txt $familySummaryFile.xls`;
#}
else
{ 
  foreach my $family(@Families)
  {
    my @outFile = split(/\//, $VARcmd{'inFile'});    
    my $familySummaryFile="$VARcmd{'outputDir'}/$family/$VARcmd{'version'}/$outFile[$#outFile].$family";    
    open(OUT, ">$familySummaryFile.txt") || die "Cant create the family summary file $familySummaryFile.txt\n";
    open(OUT_Exonic, ">$familySummaryFile.Exonic.txt") || die "Cant create the family exonic file $familySummaryFile.Exonic.txt\n";
    
    my $startCol;
    my $endCol;
    my $uniqueCol;
    my $varTagCol;
    my @familyCols;

    foreach my $line(@data)
    {           
      if($line=~/Family_ID/)
      {
        my @ss=split(/\t/, $line);
        
        for(my $i=0;$i<=$#ss;$i++)
        {
          #if($ss[$i]=~/$family.CASE/){$startCol=$i;}
          if($ss[$i]=~/_$family$/){push(@familyCols, $i)}
          if($ss[$i]=~/InAllCase/){$uniqueCol=$i}
          if($ss[$i]=~/CHROM_POS_REF_ALT/){$varTagCol=$i}
        }
        $endCol=$startCol+3;
        print OUT join("\t", @ss[0..$varTagCol]), "\t";
        map{print OUT "@ss[$_]\t"} @familyCols;
        print OUT join("\t", @ss[$uniqueCol..$#ss]), "\n";

        print OUT_Exonic join("\t", @ss[0..$varTagCol]), "\t";
        map{print OUT_Exonic "@ss[$_]\t"} @familyCols;
        print OUT_Exonic join("\t", @ss[$uniqueCol..$#ss]), "\n";
      }
      else
      {
        if($line=~/^$family\t/)
        {
          my @ss=split(/\t/, $line);
          print OUT join("\t", @ss[0..$varTagCol]), "\t";
          map{print OUT "@ss[$_]\t"} @familyCols;
          print OUT join("\t", @ss[$uniqueCol..$#ss]), "\n";

          if($line=~/\t(exonic|splicing)/ && $line!~/\t(ncRNA|synonymous|UTR|upstream)/){
            my @ss=split(/\t/, $line);
            print OUT_Exonic join("\t", @ss[0..$varTagCol]), "\t";
            map{print OUT_Exonic "@ss[$_]\t"} @familyCols;
            print OUT_Exonic join("\t", @ss[$uniqueCol..$#ss]), "\n";
          }
        }
      }
    }
    #`perl $VARcmd{'pipelineDir'}/csv2xlx.pl $familySummaryFile.txt $familySummaryFile.xls`;
  }
}
