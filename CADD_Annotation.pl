#!/usr/bin/perl 
use strict;
use Getopt::Long;

# Get command line arguments
my $outputFileSufix="CADD.txt";

my ($inFile, $pipelineDir, $baseDIR, $chr, $variant, $PromoterFile);

GetOptions ( "inFile=s"      => \$inFile,
             "pipelineDIR=s" => \$pipelineDir,
             "baseDIR=s"     => \$baseDIR,
             "chr=s"         => \$chr,
             "var=s"         => \$variant,
             "promoterBED=s" => \$PromoterFile
           );

########## Create tmp vcf file 
open(IN, "<$baseDIR/$inFile.txt") || die "cant open the file input File $inFile.txt\n";

my $chrNO;
my $posNO;
my $refNO;
my $altNO;
my $pidNO;
my $variantNO;
my $vartagNO;

my %data_GERMLINE;
my $headerPresent = 0;
my $header;
my %dataVCF;

open(TMP, ">$baseDIR/$inFile.$chr.$variant.tmp.vcf") || die "cant create the tmp file $inFile.tmp.vcf\n";

while(<IN>)
{
  chomp;
  my $line=$_;

  if($line=~/Family_ID/)
  {
    $headerPresent=1;
    $header=$line;

    my @ss=split(/\t/, $line);
    for(my $i=0;$i<=$#ss; $i++)
    {
      if($ss[$i]=~/VARIATION/){$variantNO=$i;}
      if($ss[$i]=~/CHROM/){$chrNO=$i;}
      if($ss[$i]=~/POS/){$posNO=$i;}
      if($ss[$i]=~/REF/){$refNO=$i;}
      if($ss[$i]=~/ALT/){$altNO=$i;}
    }
    print TMP "#CHROM\tPOS\tID\tREF\tALTt\tQUAL\tFILTER\tINFO\n";
  }
  else
  {
    if($headerPresent==0)
    {
      die "Header info missing in the $inFile.txt\n";
    }

    my @ss=split(/\t/, $line);
    my @VarTag = ($ss[$chrNO], $ss[$posNO], $ss[$refNO], $ss[$altNO]);
    my $key=join("_", $ss[$chrNO], $ss[$posNO], $ss[$refNO], $ss[$altNO]);

    if($VarTag[0]=~/^$chr$/ && $ss[$variantNO]=~/^$variant$/)
    {
      $dataVCF{$key}="$VarTag[0]\t$VarTag[1]\t.\t$VarTag[2]\t$VarTag[3]\t.\t.\t.";
    }
  }
}

foreach my $key (keys %dataVCF)
{
  print TMP "$dataVCF{$key}\n";
}

close TMP;
close IN;

################# 

# Annotation temp vcf with CADD
print "Running CADD ...\n";
my @CADD_data;

system("bgzip -f $baseDIR/$inFile.$chr.$variant.tmp.vcf");
system("/ibios/tbi_cluster/13.1/x86_64/CADD_v1.3/bin/score.sh $baseDIR/$inFile.$chr.$variant.tmp.vcf.gz $baseDIR/$inFile.$chr.$variant.tmp.out.vcf.gz");
chomp(@CADD_data=`zcat $baseDIR/$inFile.$chr.$variant.tmp.out.vcf.gz | grep -v "#"`);
#if($variant=~/Indel/)
#{
#  chomp(@CADD_data=`perl $pipelineDir/tools_CADD_annotator.pl $baseDIR/$inFile.$chr.$variant.tmp.vcf InDels.tsv.gz`);
#}
#else
#{
#  chomp(@CADD_data=`perl $pipelineDir/tools_CADD_annotator.pl $baseDIR/$inFile.$chr.$variant.tmp.vcf whole_genome_SNVs.tsv.gz`);
#}


#####################
# processing CADD annotations 
my %caddHash;
my %promoterHash;
my $join_cadd_Header;
my $exacsubPopAF_col;

foreach my $cadd (@CADD_data)
{
   #$cadd =~ s/\r//g;
   if($cadd=~/^#CHROM/)
   {
   }
   else
   {
     my @split_cadd = split(/\t/, $cadd);
     my $key = join("_", $split_cadd[0], $split_cadd[1], $split_cadd[2], $split_cadd[3]);
     
     $promoterHash{$key}=InPromoter($key);
     $caddHash{$key} = $split_cadd[5];
   }
}  

##############
## open Cohort combined file
print "Combinging data \n";
open(IN, "<$baseDIR/$inFile.txt") || die "cant read the ExAC annotated file $inFile.txt\n";
open(OUT, ">$baseDIR/$inFile.$chr.$variant.$outputFileSufix") || die "cant create the CADD annotation file $inFile$outputFileSufix\n";

my %data;
my $oldHeader;

while(<IN>)
{
  chomp;
  my $line = $_;
  if($line=~/Family_ID/)
  {
    print OUT "$line\tCADD_PHRED_LOCAL\tPromoterGene\n";
  }
  else
  {
    my @ss=split(/\t/, $line);
    my $key=join("_", $ss[$chrNO], $ss[$posNO], $ss[$refNO], $ss[$altNO]);
    my @VarTag=($ss[$chrNO], $ss[$posNO], $ss[$refNO], $ss[$altNO]);
    
    #### Annotation with CADD and local control
    if($VarTag[0]=~/^$chr$/ && $ss[$variantNO]=~/^$variant$/)
    {
      # -- Printing the results to ouput file 

      if(defined $caddHash{$key})
      {
        print OUT "$line\t";
        print OUT "$caddHash{$key}\t";
        print OUT "$promoterHash{$key}\n";
        
      }
      else
      {
        print OUT "$line\t";
        print OUT "NoCADD\t";
        print OUT "$promoterHash{$key}\n";
      }
    }
  }
}
close IN;

#### 
# Cleaning up
`rm $baseDIR/$inFile.$chr.$variant.tmp.vcf.gz`;
`rm $baseDIR/$inFile.$chr.$variant.tmp.out.vcf.gz`;

### Promoter 
sub InPromoter
{
  my @var=split("_", $_[0]);
  my $promoterBED=$PromoterFile;
  chomp(my @promoters = `tabix $promoterBED $var[0]:$var[1]-$var[1] | cut -f14`);
  my $genePromoter="NotInPromoter";

  if($promoters[0]!~/^$/)
  {
    $genePromoter=join(":", @promoters);
  }
  return($genePromoter);
}
