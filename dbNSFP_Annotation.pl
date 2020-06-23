#!/usr/bin/perl 
use strict;
use Getopt::Long;

my %VARcmd;

# Get command line arguments
GetOptions ( "inputFilePrefix=s"       => \$VARcmd{'inputFilePrefix'},
             "baseDIR=s" => \$VARcmd{'baseDIR'},
             "pipelineDIR=s"     => \$VARcmd{'pipelineDIR'},
             "outputFileSufix=s" => \$VARcmd{'outputFileSufix'},
             "populationBasedFilter=f" => \$VARcmd{'populationBasedFilter'},
             "ExACsubPopulation=s" => \$VARcmd{'ExACsubPopulation'},
             "MAF=f" => \$VARcmd{'MAF'},
             "chr=s"             => \$VARcmd{'chr'},
             "var=s"             => \$VARcmd{'var'}

           );

my $inFile="$VARcmd{'baseDIR'}/$VARcmd{'inputFilePrefix'}.$VARcmd{'chr'}.$VARcmd{'var'}.CADD.GTEx";
my $pipelineDir=$VARcmd{'pipelineDIR'};
my $outputFileSufix=$VARcmd{'outputFileSufix'};

########## Create tmp vcf file 
open(IN, "<$inFile.txt") || die "cant open the file input File $inFile.txt\n";

my $chrNO;
my $posNO;
my $refNO;
my $altNO;
my $pidNO;
my $annovarNO;
my $variationNO;

my %data_GERMLINE;
my $headerPresent = 0;
my $header;
my %dataVCF_exonic;
my %dataVCF_nonexonic;

open(TMP_Exonic, ">$inFile.exonic.tmp.vcf") || die "cant create the tmp file $inFile.exonic.tmp.vcf\n";
open(TMP_nonExonic, ">$inFile.nonexonic.tmp.vcf") || die "cant create the tmp file $inFile.nonexonic.tmp.vcf\n";

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
      if($ss[$i]=~/ANNOVARAnnotation/){$annovarNO=$i;}
      if($ss[$i]=~/VARIATION/){$variationNO=$i;}
      if($ss[$i]=~/PatientPIDs/){$pidNO=$i;}
    }
    print TMP_Exonic "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
  }
  else
  {
    if($headerPresent==0)
    {
      die "Header info missing in the $inFile.txt\n";
    }

    my @ss=split(/\t/, $line);
    
    if($ss[$annovarNO]=~/exonic/)
    {
      my $key = join("_", $ss[$chrNO],$ss[$posNO],$ss[$refNO],$ss[$altNO]);
      $dataVCF_exonic{$key}="$ss[$chrNO]\t$ss[$posNO]\t.\t$ss[$refNO]\t$ss[$altNO]\t.\t.\t.";
    }
    else
    {
      my $key = join("_", $ss[$chrNO],$ss[$posNO],$ss[$refNO],$ss[$altNO]);
      $dataVCF_nonexonic{$key}="$ss[$chrNO]\t$ss[$posNO]\t.\t$ss[$refNO]\t$ss[$altNO]\t.\t.\t.";
    }
  }
}

foreach my $key (keys %dataVCF_exonic)
{
  print TMP_Exonic "$dataVCF_exonic{$key}\n";
}


close TMP_Exonic;
close TMP_nonExonic;
close IN;

################
# Annotation with GTEx class 
chomp(my @gtex = `bedtools intersect -a $inFile.exonic.tmp.vcf -b GTEx/GTEx_Analysis_V4_RNA-seq_RNA-SeQCv1.1.8_exon_reads.75Percentile.80PerSamples.txt.SUMMARY.bed -wb`);

my %GTEx;

foreach my $class(@gtex)
{
  my @ss = split(/\t/, $class);
  my $key = join("_", $ss[0], $ss[1], $ss[3], $ss[4]);
  $GTEx{$key}="$ss[12]\t$ss[13]";
}


################# 

# Annotation temp vcf with dbNSFP

chomp(my @dbNSFP_data=`perl $pipelineDir/tools_dbNSFP_annotator.pl $inFile.exonic.tmp.vcf`);


#####################
# print dbNSFP annotations 
my %dbnsfpHash;
my $join_dbnsfp_Header;
my $exacsubPopAF_col;
my $dbnsfp_colNu = 0;

foreach my $dbnsfp (@dbNSFP_data)
{
   $dbnsfp =~ s/\r//g;
   if($dbnsfp=~/^#CHROM/)
   {
     my @ss=split(/\t/, $dbnsfp);
     my @dbnsfp_Header = splice(@ss, 19);
     $dbnsfp_colNu = $#dbnsfp_Header;
     for(my $i=0;$i<=$#dbnsfp_Header;$i++)
     {
       if($dbnsfp_Header[$i]=~/$VARcmd{'ExACsubPopulation'}/){$exacsubPopAF_col = $i}
     }
     $join_dbnsfp_Header = join("\t", @dbnsfp_Header);
   }
   else
   {
     my @split_dbnsfp = split(/\t/, $dbnsfp);
     my $key = join("_", $split_dbnsfp[0], $split_dbnsfp[1], $split_dbnsfp[3], $split_dbnsfp[4]);
     my @dbnsfp_Annotation = splice(@split_dbnsfp, 19);

     my $exacsubPopAF = sprintf("%.5f", $dbnsfp_Annotation[$exacsubPopAF_col]);

     ## If we need to filter for the population based data
     if($VARcmd{'populationBasedFilter'} == 1)
     {
       if($exacsubPopAF <= $VARcmd{'MAF'})
       {
         my $join_dbnsfp_Annotation = join("\t", @dbnsfp_Annotation);
         $dbnsfpHash{$key} = $join_dbnsfp_Annotation;
       }
     }
     ## otherwise no population based filter
     else
     {
       my $join_dbnsfp_Annotation = join("\t", @dbnsfp_Annotation);
       $dbnsfpHash{$key} = $join_dbnsfp_Annotation;
     }
   }
}  
     
##############
## open Cohort combined file
open(IN, "<$inFile.txt") || die "cant read the ExAC annotated file $inFile.txt\n";
open(OUT, ">$inFile.$outputFileSufix") || die "cant create the dbNSFP annotation file $inFile$outputFileSufix\n";

my %data;
my $oldHeader;

while(<IN>)
{
  chomp;
  my $line = $_;
  if($line=~/VAR_TYPE/)
  {
    print OUT "$line\tGTEx_D\tGTEx_CLASS\t$join_dbnsfp_Header\n";     
  }
  else
  {
    my @ss=split(/\t/, $line);
    my $key=join("_", @ss[7..10]);  
    # -- Printing the results to ouput file 

    # printing the basic info
    print OUT "$line\t";

    # printing the exonic variants
    if(defined $dbnsfpHash{$key})
    {
      if(defined $GTEx{$key})
      {
        print OUT "$GTEx{$key}\t$dbnsfpHash{$key}\n";
      }
      else
      {
        print OUT "NA\tNA\t$dbnsfpHash{$key}\n";
      }
    }
    # printing the non-exonic variants
    else
    {
      print OUT "NA\tNA\t";
      print OUT "\t." x $dbnsfp_colNu, "\n";
    }
  }
}
close IN;

## Cleaning up VCF file
`rm $inFile.exonic.tmp.vcf`;
`rm $inFile.nonexonic.tmp.vcf`;
