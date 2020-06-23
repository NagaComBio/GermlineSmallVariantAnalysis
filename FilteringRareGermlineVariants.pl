#!/usr/bin/env perl
#
# Program to filter the variants based on different disease models
#  Denovo variant
#  Homozygous variant
#  Compund heterzygous variants 
#  Heterozygous variants 
#  
# Can analysis TRIO and single patients ! 
#   TRIO = child, mother, father
# 
# Uses VCF file
#   SNVs calls from mpileup
#   Indel calls from Platypus
#
# if TRIO sample and the samples are with single SNVs calling,
#   use addBamInfo pipeline to add parents info to the VCF files.    
# 
# VCF with annotation information
#   dbSNPs
#   1000 genome
#   ANNOVAR
# 
# 
########################################################################
########################################################################
=a /README/

INPUT_DATA:
===========
  SNVs:  Accepts the SNV calling from platypus pipeline.  
  Indel: Accepts the Indel calling from Platypus pipeline

What_it_can_do ? 
================
  TRIO: Can take in 3 person for trio analysis
  Single: Can run on only 1 person 
  Cancer: Can compare germline and cancer samples of same patient

BeforeStaring TO_DO:
====================
  Basic_Information
  ------------------
    Create a text file with some basic information about patient
    (PID, FAMILY_ID, CASEorCONTROL are mandatory)
    Results at the end will be combined based on the Family_id. 

    For TRIO:  mother's and father's PID should be given in the same line as patient.
    For Single: same Family_id should be given to all the PID, if you want to combine the anlysis at the end. 

    Use below format (TAB separated):
     #PID FAMILY_ID CASE/CONTROL GENDER MOTHER_ID FATHER_ID PROJECT_ID
     PID_01  FAMILy_1  CASE  MALE  PID_02  PID_03  SAMPLE
     PID_02  FAMILY_1  CONTROL FEMALE  NA  NA  SAMPLE
     PID_03  FAMILY_1  CASE  MALE  NA  NA  SAMPLE
     PID_04  FAMILY_2  CASE  NA  NA  NA  SAMPLE

SOME_INFO:
==========
  Local control:
  --------------
    Removing the variants with same genotype in local control reduces false positive variants.

  BlackList:
  ----------
    Black list files are available from cancer whole genome.
  
  Intolerance Score:
  ------------------
    The intolerance score is a studentized residual score derived from all the mutations in a gene and functional mutations with MAF > 1%.

    Negative the score, the gene is highly intolerant to a new functional mutation and positive the score, it is tolerant to a new functional mutation.

  Additional_Annotation:
  ----------------------
    Known disease association of variant and the gene.
    BIOGPS
    Subcellular localization
    ExAC
    ESP_6500    

OUTPUT_DATA:
============
  Generates a list of variants filtered with different genetic models and combines them into single file for each pids.

  Filter files generated for each sample will be combined and place in the ANALYSIS_DIR for further analysis.

  A cohort summary will be generated for all the PIDs used for the particular run.
=cut
#============================================================================================
use strict;
#use warnings;
use lib '/home/paramasi/perl5/lib/perl5'; 
use Config::Simple;
use Getopt::Long;

my $GunzipError;
use IO::Uncompress::Gunzip qw($GunzipError);

# Input arguments
#  VARcmd = CmdLine arguments, PIDs and runtime specific
#  VARcnf = Stable arguments from configuration file
#  VAR    = Will contain all the input arguments

my %VARcmd;

GetOptions ( "configFile=s"      => \$VARcmd{'configFile'},
             "inputFile=s"       => \$VARcmd{'inputFile'},
             "summaryFile=s"     => \$VARcmd{'summaryFile'},             
             "patientPIDs=s"     => \$VARcmd{'patientPIDs'},             
             "variation=s"       => \$VARcmd{'variation'},
             "diseaseModel=s"    => \$VARcmd{'D_model'});


my $cnf= new Config::Simple($VARcmd{'configFile'}) or die Config::Simple->error();
my %VARcnf = $cnf->vars();

my %VAR1 = (%VARcmd, %VARcnf);
my (%VAR) = CheckInputParameters(\%VAR1);

#-----------------------------------------------------------------------
# Opening input file, creating summary and log files

my $IN;

if($VAR{'inputFile'}=~/vcf.gz$/)
{
  #$IN = IO::Uncompress::Gunzip->new( $VAR{'inputFile'} ) or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
  open($IN , 'gzip -dc '.$VAR{'inputFile'}.' |') || die "Can't open the input VCF: $VAR{'inputFile'} $!\n";
}
elsif($VAR{'inputFile'}=~/vcf$/)
{
  open ($IN, "$VAR{'inputFile'}") || die "Can't open the input VCF: $VAR{'inputFile'} $!\n";
} 
 
open (SUMMARY, ">$VAR{'summaryFile'}") || die "Cant create the output file:$VAR{'summaryFile'} $!\n";
open (LOG, ">$VAR{'summaryFile'}.log") || die "Cant craete the log file $!\n";

#-----------------------------------------------------------------------
# Simple Statistics on each step of the filtering
my $c_Allread_gt10           = 0;
my $c_SNV_more_confidence_th = 0;
my $c_SNV_match_1KG          = 0;
my $c_SNV_match_1KG_AF_1     = 0;
my $c_SNV_nt_match           = 0;
my %c_SNV_filtered           = 0;
my $c_Finalfiltered          = 0;
my $c_SNV_for_filtering      = 0;
my $c_SNV_nsy_spl            = 0;
my $c_SNV_SYS                = 0;
my $c_SNV_Non_SYS            = 0;
my $c_ExAC_common            = 0;
my $c_gnomAD_common          = 0;
my $c_LC2_common             = 0;
my $c_LC3_common             = 0;
my $c_LC4_common             = 0;
my $c_phase1_common          = 0;
my $c_phase3_common          = 0;
my $c_SNV_match_Common       = 0;
my $c_clinsig_rare           = 0;

#-----------------------------------------------------------------------
## Log information
chomp(my $date=`date`);
print LOG "Patient_pids\t$VAR{'patientPIDs'}\n";
print LOG "Input file fullpath\t$VAR{'inputFile'}\n";
print LOG "Running for TRIO (==0) or single (==1)\t$VAR{'trioORsingle'}\n";
print LOG "Date stamp\t$date\n";
print LOG "SNP Exonic\t$VAR{'patternSnpEc'}\n";
print LOG "SNP ANNOVAR\t$VAR{'patternSnpAnn'}\n";
print LOG "Indel Exonic\t$VAR{'patternIndelEc'}\n";
print LOG "Indel ANNOVAR\t$VAR{'patternIndelAnn'}\n";
print LOG "Minimum Reads\t$VAR{'MinReads'}\n";
print LOG "Minimum QUAL\t$VAR{'QUAL'}\n";
print LOG "Maximum MAF\t$VAR{'RareMAF'}\n";
print LOG "Minimum VAF for Denovo variant\t$VAR{'MinVAF_HT'}\n";
print LOG "Maximum VAF for Heterozygous variant\t$VAR{'MaxVAF_HT'}\n";
print LOG "Maximum VAF for Homozygous Reference\t$VAR{'MaxVAF_HoRef'}\n";
print LOG "Maximum VAF for Homozygous Alternate\t$VAR{'MaxVAF_HoAlt'}\n";
print LOG "Genotype quality (GT)\t$VAR{'genotypeQuality'}\n\n";

print LOG "Baisc-Info\t$VAR{'basicInfo'}\n";
print LOG "GeneticAssocationFile\t$VAR{'GADdb'}\n";
print LOG "GeneDescription\t$VAR{'geneDesc'}\n";
print LOG "OMIM file\t$VAR{'omimGeneMap'}\n";
print LOG "GO annotation file\t$VAR{'GOannotation'}\n";
print LOG "SNP-Indel parentVCF file\t$VAR{'RaKi_UGD_parentVCF'}\n\t$VAR{'HIPO_EURO_parentVCF'}\n";
print LOG "Intolerance score Local\t$VAR{'RVIS_LocalControlFile'}\n";
print LOG "Intolerance score PLOS\t$VAR{'RVIS_PlosGeneticsFile'}\n";
print LOG "Intolerance score ExAC\t$VAR{'RVIS_ExAC'}\n";
print LOG "BioGPS Annotation file\t$VAR{'bioGPS_AnnFile'}\n";
print LOG "ESP_ControlFile\t$VAR{'ESP_File'}\n";
print LOG "SCL information\t$VAR{'SCLinfoFile'}\n";


# ANNOVAR exonic/functional classification patterns
my @pattern_FN          = ();  
my @query_FN            = ();
my $query_EC            = "";
my $query_AN_FN         = "";

my ($query_dbsnp, $query_1kG);
my ($motherPid, $fatherPid);

my ($child_DP, $mother_DP, $father_DP, $DP_AllAbove10or20, $control_DP, $tumor_DP);
my ($child_BaseQ, $mother_BaseQ, $father_BaseQ, @BaseQ, $control_BaseQ, $tumor_BaseQ);
my (@FowreadRevrread, $control_FowreadRevrread, $tumor_FowreadRevrread); 
my ($child_DP_field, $mother_DP_field, $father_DP_field, @DP_field, $control_DP_field, $tumor_DP_field);
my ($child_fr, $child_rr, $child_fnr,$child_rnr);
my ($mother_fr, $mother_rr, $mother_fnr,$mother_rnr);
my ($father_fr, $father_rr, $father_fnr,$father_rnr);
my ($control_fr, $control_rr, $control_fnr, $control_rnr);
my ($tumor_fr, $tumor_rr, $tumor_fnr, $tumor_rnr);
my ($child_alt_AF, $mother_alt_AF, $father_alt_AF, $control_alt_AF, $tumor_alt_AF);
my (@var_addbaminfo, $mother_addbaminfo, $father_addbaminfo);
my @query_alt_AF;
my @query_dbSNP_1KG;

# Minimum and 200% minimum reads 
my $miniReads = $VAR{'MinReads'};
my $mini200Reads = $VAR{'MinReads'} * 2;

# SNVs or Indel
my $st_variation = "";
#if($VAR{'variation'}    =~ /SNVs/){$st_variation="SNVs"}
#elsif($VAR{'variation'} =~ /Indel/){$st_variation="Indel"}

########################################################################
#
# Storing information for later use in the program
#
########################################################################

# Basic Information 
#  File where we store all the basic information about the samples
#  It should ba TAB limited txt file
#  Leave the columns as NA if the information is not available
#  PID	FAMILY_ID	CASE_CONTROL	GENDER	MOTHER_ID	FATHER_ID	PROJECT

open(BI, "<$VAR{'basicInfo'}") || die "can't open the $VAR{'basicInfo'} file\n";
my %patientGroup;
my %CaseControl;
my %gender;
my %mother;
my %father;

while(<BI>)
{
  chomp;
  my @split_BI = split(/\t/, $_);
  $patientGroup{$split_BI[0]}=$split_BI[1];
  $CaseControl{$split_BI[0]}=$split_BI[2];
  $gender{$split_BI[0]}=$split_BI[3];
  $mother{$split_BI[0]}=$split_BI[4];
  $father{$split_BI[0]}=$split_BI[5];
}

# If it is a TRIO run, mother and father id is missing in the
#  basic information file, the program dies

if(($mother{$VAR{'patientPIDs'}}=~/NA/ || $father{$VAR{'patientPIDs'}}=~/NA/) && $VAR{'trioORsingle'}==0)
{
  die "TRIO pipeline; Must be parents IDs; Or parents IDs are not available";
}
else
{
  $motherPid=$mother{$VAR{'patientPIDs'}};
  $fatherPid=$father{$VAR{'patientPIDs'}};
}
close BI;

print LOG "MotherPid\t$motherPid\n";
print LOG "FatherPid\t$fatherPid\n";
# Annotating with genetic association ##################################

open(GAD, "<$VAR{'GADdb'}") || die "Can't open the $VAR{'GADdb'} file\n";
my %GADInfo = ();
while(<GAD>)
{
  chomp;
  my $line = $_;
  my @split_GAD = split(/\t/, $_);  
  my $gene = $split_GAD[8];
  push(@{$GADInfo{$gene}},$split_GAD[2]);  
}
close GAD;

# Gene Description #####################################################

open(GD, "<$VAR{'geneDesc'}") || die "Cant open the $VAR{'geneDesc'} file\n";
my %GeneDesc = ();
while(<GD>)
{
  chomp;
  my @split_line=split(/\t/, $_);
  my @previousSymbol = map {split(/, /, $_)} $split_line[4];
  my @synonyms = map {split(/, /, $_)} $split_line[5];
  
  $GeneDesc{$split_line[1]}=$split_line[2];
  #map{$GeneDesc{$_}=$split_line[2]} @previousSymbol;
  #map{$GeneDesc{$_}=$split_line[2]} @synonyms;  
}
close GD;

# OMIM #################################################################

open(OMIM, "<$VAR{'omimGeneMap'}") || die "Cant open the $VAR{'omimGeneMap'}\n";
my %OMIM_Disorder = ();
my %OMIM_Number = ();

while(<OMIM>)
{
  chomp;
  my @split_omim = split(/\|/, $_);
  my @split_gene = split(/, /, $split_omim[5]);
  #print join("# ", @split_gene), "\t $split_gene[11]\n";
  map{$OMIM_Number{$_}=$split_omim[8]} @split_gene;
  map{$OMIM_Disorder{$_}=$split_omim[11]} @split_gene;
}
close OMIM;

# Gene Ontology ########################################################

open(GO, "<$VAR{'GOannotation'}") || die "Cant open the $VAR{'GOannotation'}\n";
my %GOinfo = ();
while(<GO>)
{
 chomp;
 my @split_GO = split(/\t/, $_);
 my @split_geneSymbol = split(/\|/, $split_GO[3]);
 my @geneSym = map{/:(.*)$/ ? $1 : ()} @split_geneSymbol;
 map{push(@{$GOinfo{$_}}, $split_GO[1])} @geneSym;
}
close GO;

# Parent Control #######################################################
#  Local exome control database
#  Generate custom localControl from samples with different phenotype

my %P_RAKI_UGD;
my %P_HIPO_EURO;

if($VAR{'RUN_LOCAL_CONTROL'} == 1)
{
  open(PRU, "<$VAR{'RaKi_UGD_parentVCF'}") || die "Cant open the $VAR{'RaKi_UGD_parentVCF'}\n";
  open(PHE, "<$VAR{'HIPO_EURO_parentVCF'}") || die "Cant open the $VAR{'HIPO_EURO_parentVCF'}\n";

  while(<PRU>)
  {
    chomp;
    my @split_PRU = split(/\t/, $_);
    my $value = pop(@split_PRU);
    my $key = join("_", @split_PRU);
    $P_RAKI_UGD{$key}=$value;
  }

  while(<PHE>)
  {
    chomp;
    my @split_PHE = split(/\t/, $_);
    my @key_1 = splice(@split_PHE, 0,4);
    my $value = shift(@split_PHE);
    my $key = join("_", @key_1);
    $P_HIPO_EURO{$key}=$value;
  }
  close PRU;
  close PHE;
}

# BlackList ############################################################
#  BlackList from Medulloblastoma WGS
#  Available

my %BlackList;

if($VAR{'RUN_BLACK_LIST'} == 1)
{  
  open(BL, "<$VAR{'BLACK_LIST_PRUNED'}") || die "Cant open the $VAR{'BLACK_LIST_PRUNED'}\n";

  while(<BL>)
  {
    chomp;
    my @split_BL = split(/\t/, $_);
    my $value = pop(@split_BL);
    my $key = join("_", @split_BL);
    $BlackList{$key}=$value;
  }   
  close BL;
}

# Intolerate scores ####################################################
#  Student residual score
#   Local exome control
#   ExAC 63,000 exome data
#   ESP 6500 exome data

open(RVIS_LC, "<$VAR{'RVIS_LocalControlFile'}") || die "Cant open the file $VAR{'RVIS_LocalControlFile'}\n";
my %RVIS_LC = ();
my %LocalControlMutation = ();
while(<RVIS_LC>)
{
  chomp;
  my @LS=split("\t", $_);
  $RVIS_LC{$LS[1]}=$LS[5];
  my $value=join("\t", $LS[2], $LS[3]);
  $LocalControlMutation{$LS[1]}=$value;
}
close RVIS_LC;

#----------------------------------------------------------------------
open(RVIS_ExAC, "<$VAR{'RVIS_ExAC'}") || die "Cant open the file $VAR{'RVIS_ExAC'}\n";
my %RVIS_ExAC = ();

while(<RVIS_ExAC>)
{
  chomp;
  my @LS=split("\t", $_);
  $RVIS_ExAC{$LS[1]}=$LS[5];
}
close RVIS_ExAC;

#-----------------------------------------------------------------------
open(RVIS_PG, "<$VAR{'RVIS_PlosGeneticsFile'}") || die "Cant open the file $VAR{'RVIS_PlosGeneticsFile'}\n";
my %RVIS_PG = ();
while(<RVIS_PG>)
{
  chomp;
  my @PG=split("\t", $_);
  $RVIS_PG{$PG[0]}=$PG[1];
}
close RVIS_PG;

# B i o G P S ##########################################################
#   Microarray expression information for top 5 tissues

open(BioGPS, "<$VAR{'bioGPS_AnnFile'}") || die "Can't open the $VAR{'bioGPS_AnnFile'} file\n";
my %bioGPS = ();
while(<BioGPS>)
{
  chomp;
  my @ss_bioGPS=split(/\t/, $_);
  $bioGPS{$ss_bioGPS[0]}=$ss_bioGPS[2];
}
close BioGPS;

# ExAC Processed Data ##################################################
my %ExACdata;
my %gnomADdata;

# S C L and Protein T Y P E ############################################
#  Data downloaded from Ingenuity knowledge base

open(SCL, "<$VAR{'SCLinfoFile'}") || die;
my %SCLinfo;
my %typeInfo;

while(<SCL>)
{
  chomp;
  my @split_SCL = split(/\t/, $_);
  $SCLinfo{$split_SCL[0]}=$split_SCL[4];
  $typeInfo{$split_SCL[0]}=$split_SCL[5];
}
close SCL;

######## ClinVar #######################################################
#open(ClinVAR, "<$VAR{'clinVar'}") || die;
my %clinVarInfo;

############ DisGeNET ##################################################
open(DGNET, "<$VAR{'DisGeNET'}") || die;
my %DisGeNETInfo;
while(<DGNET>)
{
  chomp;
  my @split_DGNET = split(/\t/, $_);
  push(@{$DisGeNETInfo{$split_DGNET[0]}}, $split_DGNET[1]);
}
close DGNET;

########################################################################
my %LocalControl_2;
my %Phase3;
my %LocalControl_3;
my %LocalControl_4;

########################################################################
#
# Let us open the input VCF file
#
########################################################################

#my @splitInputLine;

my ($patientVcfCol, $motherVcfCol, $fatherVcfCol, $controlVcfCol, $tumorVcfCol);
my ($EC_VcfCol, $AN_VcfCol, $dbsnp_VcfCol, $ThoGe_VcfCol);
my ($gene_VcfCol, $Transcript_VcfCol, $QUAL_VcfCol, $GQ_VcfCol, $GQ_Format_VcfCol);
my ($confidence_VcfCol, $classification_VcfCol);
my ($child_FowreadRevrread, $mother_FowreadRevrread, $father_FowreadRevrread);
my ($phase3_VcfCol, $ExAC_VcfCol, $gnomAD_VcfCol, $LC2_VcfCol, $LC3_VcfCol, $LC4_VcfCol, $clinvar_vcfCol);

my $headerPresent=0;

while(<$IN>)
{
  chomp;
  # Printing the comment sections to the output file
  if($_=~/^#/)
  {
    print SUMMARY "$_\n", if($VAR{'toPrintSummaryORVCF'} == 1);
      
    if($_=~/^#CHROM/)
    {
      $headerPresent=1;
      print SUMMARY "$_\n", if($VAR{'toPrintSummaryORVCF'} == 1);
      
      my @split_header=split(/\t/, $_);
      # Selecting the column number for the different sample in TRIO
      for my $i (0 .. $#split_header)
      {
        if($VAR{'trioORsingle'} == 0 )
        {
          if($VAR{'variation'} =~ /Indel|SNVs/)
          {
            $motherVcfCol  = $i, if($split_header[$i]=~/$VAR{'GERMLINE_PREFIX'}$motherPid$/i);
            $fatherVcfCol  = $i, if($split_header[$i]=~/$VAR{'GERMLINE_PREFIX'}$fatherPid$/i);
            $patientVcfCol = $i, if($split_header[$i]=~/$VAR{'GERMLINE_PREFIX'}$VAR{'patientPIDs'}$/i);
            $GQ_VcfCol = $patientVcfCol;
            #$patientVcfCol = $i, if($split_header[$i]=~/sample_/); 
            #$GQ_VcfCol = $i, if($split_header[$i]=~/sample_control/);
          }
          #elsif($VAR{'variation'} =~ /SNVs/)
          #{
          #  $motherVcfCol  = $i, if($split_header[$i]=~/^INFO_MOTHER/);
          #  $fatherVcfCol  = $i, if($split_header[$i]=~/^INFO_FATHER/);
          #  $patientVcfCol = $i, if($split_header[$i]=~/^INFO$/);          
          #  #$GQ_VcfCol         = $i, if($split_header[$i]=~/GT/);
          #}
        }
        elsif($VAR{'trioORsingle'} == 1)
        {
          if($VAR{'variation'} =~ /Indel|SNVs/)
          {
            $patientVcfCol = $i, if($split_header[$i]=~/$VAR{'GERMLINE_PREFIX'}$VAR{'patientPIDs'}$VAR{'GERMLINE_SUFFIX'}$/i);
            $GQ_VcfCol = $patientVcfCol;
          }
        }
        # Special case of cancer sample
        elsif($VAR{'trioORsingle'} == 2)
        {
	  if($VAR{'variation'} =~ /SNVs|Indel/ && $VAR{'D_model'} !~/Somatic/)
	  {
	    $controlVcfCol = $i, if($split_header[$i]=~/$VAR{'GERMLINE_PREFIX'}$VAR{'patientPIDs'}$VAR{'GERMLINE_SUFFIX'}$/i);
	    $tumorVcfCol   = $i, if($split_header[$i]=~/$VAR{'SOMATIC_PREFIX'}$VAR{'patientPIDs'}$VAR{'SOMATIC_SUFFIX'}$/i);
            $GQ_VcfCol = $controlVcfCol; 
	  }          
          elsif($VAR{'variation'} =~ /SNVs|Indel/ && $VAR{'D_model'}=~/Somatic/)
          {
            $controlVcfCol = $i, if($split_header[$i]=~/^INFO_control/);
            $tumorVcfCol   = $i, if($split_header[$i]=~/^INFO$/);
            $GQ_VcfCol = $i, if($split_header[$i]=~/$VAR{'patientPIDs'}|tumor_/);
          }
          #else
          #{
          #  $controlVcfCol = $i, if($split_header[$i]=~/blood_|control_|buffy_coat/);
          #  $tumorVcfCol   = $i, if($split_header[$i]=~/tumor_|metastasis_/);
          #  $GQ_VcfCol = $controlVcfCol;
          #}
	}
        # General columns
	$GQ_Format_VcfCol  = $i, if($split_header[$i]=~/^FORMAT$/);
        $QUAL_VcfCol       = $i, if($split_header[$i]=~/^QUAL$/);
        $EC_VcfCol         = $i, if($split_header[$i]=~/EXONIC_CLASSIFICATION/);
        $AN_VcfCol         = $i, if($split_header[$i]=~/ANNOVAR_FUNCTION/);
        $dbsnp_VcfCol      = $i, if($split_header[$i]=~/DBSNP/);
        $ThoGe_VcfCol      = $i, if($split_header[$i]=~/^1K_GENOMES$/);
        $gene_VcfCol       = $i, if($split_header[$i]=~/^GENE$/);
        $Transcript_VcfCol = $i, if($split_header[$i]=~/ANNOVAR_TRANSCRIPTS/);
        $confidence_VcfCol = $i, if($split_header[$i]=~/^CONFIDENCE$/);
        $classification_VcfCol = $i, if($split_header[$i]=~/CLASSIFICATION$/);
        $phase3_VcfCol     = $i, if($split_header[$i]=~/PHASE3_1K_GENOMES/);
        $ExAC_VcfCol       = $i, if($split_header[$i]=~/gnomAD_EXOMES_v2.1/);
	$gnomAD_VcfCol     = $i, if($split_header[$i]=~/gnomAD_GENOMES_v2.1/);
        $LC2_VcfCol     = $i, if($split_header[$i]=~/LOCAL_CONTROL_2/);
        $LC3_VcfCol     = $i, if($split_header[$i]=~/LOCAL_CONTROL_3/);
        $LC4_VcfCol     = $i, if($split_header[$i]=~/LOCAL_CONTROL_4/);
        $clinvar_vcfCol = $i, if($split_header[$i]=~/Clinvar/);
      }
      # Checking for the column information   
      if($GQ_VcfCol=~/^$/){die "Genotype info is missing in VCF file \n";}
      if($QUAL_VcfCol=~/^$/){die "Base qual column is missing in the VCF file\n";}
      if($EC_VcfCol =~/^$/){die "Exonic classification is missing in VCF file\n";}
      if($AN_VcfCol =~/^$/){die "ANNOVAR column is missing in VCF file\n";}
      if($dbsnp_VcfCol =~/^$/){die "dbSNP annotation is missing in VCF file\n";}
      if($ThoGe_VcfCol =~/^$/){die "1000 genome annotation is missing in VCF file\n";}
      if($gene_VcfCol =~/^$/){die "Gene Symbol/Name is missing in VCF file\n";}
      if($Transcript_VcfCol =~/^$/){print "Transcipt information is missing in VCF file\n";}
      
      #Check for parent column, if its for TRIO run
      if($VAR{'trioORsingle'} == 0)
      {
        if($fatherVcfCol =~/^$/){die "Father genotype is missing in VCF file\n";}
        if($motherVcfCol =~/^$/){die "Mother genotype is missing in VCF file\n";}  
      }

      # Checking for Somatic Model 
      if($VAR{'trioORsingle'}==2)
      {
        if($controlVcfCol=~/^$/){die "No control info data\n"};
	if($tumorVcfCol=~/^$/){die "No tumor info data\n"};
      
        if($VAR{'D_model'}=~/Somatic/) # If the confidenceScore and Classification based filtering is on
        {
          if($confidence_VcfCol=~/^$/){die "Confidence annotation is missing in the VCF file\n";}
          if($classification_VcfCol=~/^$/){die "Variant classification is missing in the VCF file\n";}
        }        
      }
      else
      {
        if($patientVcfCol =~/^$/){die "Patient INFO or Genotype column missing in VCF file\n";}
      }
      
      if($VAR{'toPrintSummaryORVCF'} == 0)	  
      {
        print SUMMARY "VAR_TYPE\tFamily_ID\tPatientPIDs\tORIGIN\tGeneticModel\tVARIATION\tdbSNPsID\t";
        print SUMMARY "CHROM\tPOS\tREF\tALT\t";
        print SUMMARY "GENE\tDP_AllAbove10or20\t";
        if($VAR{'trioORsingle'}==2)
        {
	  print SUMMARY "ControlINFO\tTumorINFO\tEMPTYINFO\t";	
          #print SUMMARY "ControlBaseQ\tTumorBaseQ\tEMPTYBaseQ\t";
          print SUMMARY "ControlVAF\tTumorVAF\tEMPTYVAF\t";
          print SUMMARY "ControlGT\tTumorGT\tEMPTYGT\t";
        }
        else
        {
          print SUMMARY "ChildINFO\tMotherINFO\tFatherINFO\t";
	  #print SUMMARY "ChildBaseQ\tMotherBaseQ\tFatherBaseQ\t";
	  print SUMMARY "ChildVAF\tMotherVAF\tFatherVAF\t"; 
	  print SUMMARY "ChildGT\tMotherGT\tFatherGT\t";
        }
        print SUMMARY "EXONIC_CLASSIFICATION\tANNOVARAnnotation\t";
        print SUMMARY "TRANSCRIPT\tdbSNPAnnotation\t1000GenomeAnnotation\t";
        print SUMMARY "gnomAD_Exomes_MAF\tgnomAD_Exomes_Homo\tgnomAD_Exomes\t";
        print SUMMARY "gnomAD_Genomes_MAF\tgnomAD_Genomes_Homo\tgnomAD_Genomes\t";
        print SUMMARY "LocalControl2_AF\t";
        print SUMMARY "LocalControl3_AF\t1kG_Phase3_AF\t1KG_Phase3\tLocalControl4_AF\t";
        print SUMMARY "Rareness\tClinvarness\t";
        #print SUMMARY "inBlackList\tBlackList\t";
        #print SUMMARY "inLocalControl_1\tLocalControl_1\t";
        #print SUMMARY "inLocalControl_2\tLocalControl_2\t";
        #print SUMMARY "InAll_LocalControl-1-2-BlackList\t";
        print SUMMARY "GeneticAssociationDB\tDiseaseGeneDB\t";
        print SUMMARY "OMIM_ID\tOMIM_Disorder\tClinVAR\tClinVAR_significance\t";
        print SUMMARY "AllVariants\tAllFunctionalVariantsWithMAF_>_1%\t";
        print SUMMARY "IntoleranceScore_LocalControl\tIntoleranceScore_PlosGenetics\t";
        print SUMMARY "IntoleranceScore_ExAC\t";
        print SUMMARY "GeneDesc\tGO_Biological_Process\tBioGPS\tSCL\tProteinType\n";
        #print SUMMARY "child/Control_ReadInFowardReverse\tmother/tumor_ReadInFowardReverse\tfather_ReadInFowardReverse\n";
        #print SUMMARY "CADD_PHRED_LOCAL\t";
        #print SUMMARY "GENE_Promoter\n";
      }  
    }
  }
  # Loop through each variant in the VCF file
  else 
  {    
    if($headerPresent==0){die "Header is missing in $VAR{'inputFile'} \n"}  
    my @splitInputLine=split(/\t/, $_);
    my $varTag;
    
    # Variants length
    if(length($splitInputLine[3]) ==1 && length($splitInputLine[4]) == 1)
    {
      $VAR{'variation'} ="SNVs";
      $st_variation="SNVs";
    }
    else
    {
      $VAR{'variation'} ="Indel";
      $st_variation="Indel";
    }

    ### Multi or signle allele in alternate, talking only first allele, in case of multi
    if($splitInputLine[4]=~/,/)
    {
      my @alt = split(/,/, $splitInputLine[4]);
      $varTag = join("_",  @splitInputLine[0..1], $splitInputLine[3], $alt[0]);
    }
    else
    {
      $varTag = join("_", @splitInputLine[0..1], @splitInputLine[3..4]);
    }
   
    
    $query_dbsnp = $splitInputLine[$dbsnp_VcfCol];
    $query_1kG   = $splitInputLine[$ThoGe_VcfCol];
      
    $query_EC    = $splitInputLine[$EC_VcfCol];
    $query_AN_FN = $splitInputLine[$AN_VcfCol];
    
    @query_dbSNP_1KG = ($query_dbsnp, $query_1kG);
    @query_FN        = ($query_EC, $query_AN_FN);
      
    if($VAR{'trioORsingle'} == 0)
    {    
      $child_DP_field  = $splitInputLine[$patientVcfCol];
      $father_DP_field = $splitInputLine[$fatherVcfCol];
      $mother_DP_field = $splitInputLine[$motherVcfCol];
      
      if($VAR{'variation'} =~/SNVs/) # SNVs
      {
        @pattern_FN = ($VAR{'patternSnpEc'}, $VAR{'patternSnpAnn'});      
      }
      elsif($VAR{'variation'} =~/Indel/) # Indel
      {
        @pattern_FN = ($VAR{'patternIndelEc'}, $VAR{'patternIndelAnn'});        
      }      
      
      # Sent to subroutine to parse the DP field
      #  check subroutine for explanation on arguments

      ($child_DP, $child_alt_AF, $child_BaseQ, $child_FowreadRevrread)    = Extract_DP($child_DP_field, $VAR{'variation'}, $VAR{'D_model'}, $VAR{'trioORsingle'});
      ($mother_DP, $mother_alt_AF, $mother_BaseQ, $mother_FowreadRevrread) = Extract_DP($mother_DP_field, $VAR{'variation'}, $VAR{'D_model'}, $VAR{'trioORsingle'});
      ($father_DP, $father_alt_AF, $father_BaseQ, $father_FowreadRevrread) = Extract_DP($father_DP_field, $VAR{'variation'}, $VAR{'D_model'}, $VAR{'trioORsingle'});
     
      if($child_DP >= $mini200Reads && $mother_DP >= $mini200Reads && $father_DP >= $mini200Reads)
      {
        $DP_AllAbove10or20 = $mini200Reads;
      }
      elsif($child_DP >= $miniReads && $mother_DP >= $miniReads && $father_DP >= $miniReads )  
      {
        $DP_AllAbove10or20 = $miniReads;
      } 
            
      @DP_field     = ($child_DP_field, $mother_DP_field, $father_DP_field);
      @query_alt_AF = ($child_alt_AF, $mother_alt_AF, $father_alt_AF);
      @BaseQ        = ($child_BaseQ, $mother_BaseQ, $father_BaseQ);
      @FowreadRevrread =($child_FowreadRevrread, $mother_FowreadRevrread, $father_FowreadRevrread);
    }
    # Single call
    # Single patient
    elsif($VAR{'trioORsingle'} == 1)
    { 
      $child_DP_field = $splitInputLine[$patientVcfCol];
  
      if($VAR{'variation'} =~/SNVs/) # SNVs
      {      
        @pattern_FN = ($VAR{'patternSnpEc'}, $VAR{'patternSnpAnn'});  
      }
      elsif($VAR{'variation'}=~/Indel/) # Indels
      {
        @pattern_FN = ($VAR{'patternIndelEc'}, $VAR{'patternIndelAnn'});    
      }         
            
      ($child_DP, $child_alt_AF, $child_BaseQ, $child_FowreadRevrread) = Extract_DP($child_DP_field, $VAR{'variation'}, $VAR{'D_model'}, $VAR{'trioORsingle'});

      if($child_DP >= $mini200Reads){$DP_AllAbove10or20 = $mini200Reads}
      elsif($child_DP >= $miniReads){$DP_AllAbove10or20 = $miniReads;}
      
      
      @DP_field     = ($child_DP_field, 'NA', 'NA');
      @query_alt_AF = ($child_alt_AF, 'NA', 'NA');
      @BaseQ        = ($child_BaseQ, 'NA', 'NA');
      @FowreadRevrread =($child_FowreadRevrread, 'NA', 'NA');
    }
    # Cancer Samples    
    #  Control data is considered as 'Child'
    #   Tumor data is considered as 'Mother' 
    elsif($VAR{'trioORsingle'} == 2)
    {
      $control_DP_field = $splitInputLine[$controlVcfCol];
      $tumor_DP_field = $splitInputLine[$tumorVcfCol];

      if($VAR{'variation'} =~/SNVs/) # SNVs
      {      
        @pattern_FN = ($VAR{'patternSnpEc'}, $VAR{'patternSnpAnn'});  
      }
      elsif($VAR{'variation'}=~/Indel/) # Indels
      {
        @pattern_FN = ($VAR{'patternIndelEc'}, $VAR{'patternIndelAnn'});    
      }         
            
      ($control_DP, $control_alt_AF, $control_BaseQ, $control_FowreadRevrread) = Extract_DP($control_DP_field, $VAR{'variation'}, $VAR{'D_model'}, $VAR{'trioORsingle'});
      ($tumor_DP, $tumor_alt_AF, $tumor_BaseQ, $tumor_FowreadRevrread) = Extract_DP($tumor_DP_field, $VAR{'variation'}, $VAR{'D_model'}, $VAR{'trioORsingle'});

      if($control_DP >= $mini200Reads){$DP_AllAbove10or20 = $mini200Reads}
      elsif($control_DP >= $miniReads){$DP_AllAbove10or20 = $miniReads;}
      
      
      @DP_field     = ($control_DP_field, $tumor_DP_field, 'NA');
      @query_alt_AF = ($control_alt_AF, $tumor_alt_AF, 'NA');
      @BaseQ        = ($control_BaseQ, $tumor_BaseQ, 'NA');
      @FowreadRevrread =($control_FowreadRevrread, $tumor_FowreadRevrread, 'NA');
    }
    ####################################################################
    #
    # Send to filter
    #
    ####################################################################
    my $dbsnp = ""; 
  
    #print "$varTag\n";
    # Min DP reads for patient in the case of TRIO/Single or control for germline-cancer
    # Cancer somatic snvs or indels are allowed without DP checks    
    if($child_DP >= $VAR{'MinReads'} || ($VAR{'trioORsingle'} == 2 && $VAR{'D_model'}=~/somatic/i) || $control_DP >= $VAR{'MinReads'})
    {
      # Mother's and Father's DP are checked for TRIO samples	 
      if(($mother_DP >= $VAR{'MinReads'} && $father_DP >= $VAR{'MinReads'} && $VAR{'trioORsingle'} == 0) || $VAR{'trioORsingle'} == 1 || $VAR{'trioORsingle'} == 2 )
      {      
        $c_Allread_gt10++;
        # dbSNP 
        ($dbsnp) = $query_dbsnp=~/MATCH=exact;/ ? $query_dbsnp =~/;ID=(rs\d+)$/ : "NodbSNP";
         
        my @FilteringParameters = ($VAR{'D_model'}, \@query_FN, \@pattern_FN, $VAR{'variation'}, $VAR{'patientPIDs'}, \@DP_field, $VAR{'D_model'}, $st_variation, \@query_alt_AF, \@BaseQ, $dbsnp, \@query_dbSNP_1KG, \@splitInputLine, $splitInputLine[$GQ_VcfCol], $DP_AllAbove10or20, $splitInputLine[$GQ_Format_VcfCol], \@FowreadRevrread);
      
        # QUAL are checked only for TRIO and Single 
        # NO checks for Somatic snvs or Indels
        # Platypus filter has to be PASS for all germline-cancer 
        my $Filter_PASS;

        if($VAR{'trioORsingle'} == 2 && $VAR{'D_model'}=~/somatic/i)  {
          $Filter_PASS = ".";
        }
        elsif($VAR{'trioORsingle'} == 2) {
          $Filter_PASS = "(PASS|alleleBias|QD;alleleBias)";
        }
        else {
          $Filter_PASS = "PASS";
        }


        my ($hapScore)=$splitInputLine[7]=~/HapScore=(\d+);MGOF/;
        if($hapScore <= 8 && $splitInputLine[6]=~/^HapScore$/){$splitInputLine[6]="PASS"}

        if(($splitInputLine[$QUAL_VcfCol] > $VAR{'QUAL'} && $splitInputLine[6] =~/^$Filter_PASS$/ && ($VAR{'trioORsingle'} == 0 || $VAR{'trioORsingle'} == 1 || ($VAR{'trioORsingle'} && $VAR{'D_model'}=~/Germline/))) || ($VAR{'trioORsingle'} == 2 && $VAR{'D_model'}=~/somatic/i))
        {
          $c_SNV_more_confidence_th++;                            
          
          # Annotation of Rareness is checked from ExAC, 1000 genomes and dbSNP 
          # ESP-6500 is not used, becoz the issue of major-minor allele  


          # Tabix the files
          my ($exac_AF, $exac_homo, $exac_allpop) = ExAC_sub($splitInputLine[$ExAC_VcfCol], $VAR{'RareMAF'});
          my ($gnomAD_AF, $gnomAD_homo, $gnomAD_allpop) = ExAC_sub($splitInputLine[$gnomAD_VcfCol], $VAR{'RareMAF'});
          my $lc_2_af = LC_sub($splitInputLine[$LC2_VcfCol]);
          my $lc_3_af = LC_sub($splitInputLine[$LC3_VcfCol]);
          my $lc_4_af = LC_sub($splitInputLine[$LC4_VcfCol]);
          my ($Phase3_af, $phase3_allpop) =Phase3_sub($splitInputLine[$phase3_VcfCol], $VAR{'cutoff_Phase3'});
          my ($clinVar, $clinsig) = clinSig($splitInputLine[$clinvar_vcfCol]);
      
          $ExACdata{$varTag}{'MAF'}      = $exac_AF;
          $ExACdata{$varTag}{'homo'}     = $exac_homo;
          $ExACdata{$varTag}{'allpop'}   = $exac_allpop;
          $ExACdata{$varTag}{'info'}     = $splitInputLine[$ExAC_VcfCol];
          $gnomADdata{$varTag}{'MAF'}    = $gnomAD_AF;
          $gnomADdata{$varTag}{'homo'}   = $gnomAD_homo;
          $gnomADdata{$varTag}{'allpop'} = $gnomAD_allpop;
          $gnomADdata{$varTag}{'info'}   = $splitInputLine[$gnomAD_VcfCol];

          $LocalControl_2{$varTag}{'AF'} = $lc_2_af;
          $LocalControl_3{$varTag}       = $lc_3_af;
          $LocalControl_4{$varTag}       = $lc_4_af;
          $Phase3{$varTag}               = $Phase3_af;
          $clinVarInfo{$varTag}          = "$splitInputLine[$clinvar_vcfCol]\t$clinsig";

          my $Rareness="NA";
          my $clinvarness=".";
          my $seenINdb=0;
          
          # ExAC or gnomAD exomes
          if(defined $ExACdata{$varTag}{'MAF'})
          {
            #if($ExACdata{$varTag}{'MAF'} > $VAR{'RareMAF'}){$Rareness="Common"; $c_ExAC_common++;}
            if($ExACdata{$varTag}{'allpop'}=~/Common/){$Rareness="Common"; $c_ExAC_common++;}   
            $seenINdb++;
          }

          # gnomAD genomes
          if(defined $gnomADdata{$varTag}{'MAF'})
          {
            if($gnomADdata{$varTag}{'allpop'}=~/Common/){$Rareness="Common"; $c_gnomAD_common++;}
            $seenINdb++;
          }

	  ## Local control 2 - UniExome
	  if(defined $LocalControl_2{$varTag}{'AF'}) 
	  {
	    if($LocalControl_2{$varTag}{'AF'} > $VAR{'cutoff_LocalControl_2'}) {$Rareness="Common"; $c_LC2_common++;}
	    $seenINdb++;
	  }

	  ## Local control 3 - MMML
	  if(defined $LocalControl_3{$varTag})
	  {
	    if($LocalControl_3{$varTag} > $VAR{'cutoff_LocalControl_3'}){$Rareness="Common"; $c_LC3_common++;}
	    $seenINdb++;
	  }

	  ## Local control 4 - X10 samples
	  if(defined $LocalControl_4{$varTag})
	  {
	    if($LocalControl_4{$varTag} > $VAR{'cutoff_LocalControl_4'}){$Rareness="Common"; $c_LC4_common++;}
	    $seenINdb++;
	  }

	  ### Phase 1 1000 genome
	  if($query_1kG=~/MATCH/)
	  {
	    my ($allele_freq) = $query_1kG=~/;AF=(\d\.\d+);/;
	    if($allele_freq > $VAR{'cutoff_Phase1'}){$Rareness="Common"; $c_phase1_common++;}
	    $seenINdb++;
	  }           

	  ## Phase 3 1000 genome
	  if(defined $Phase3{$varTag})
	  {
	    #if($Phase3{$varTag} > $VAR{'cutoff_Phase3'}){$Rareness="Common"; $c_phase3_common++;}
            if($phase3_allpop=~/Common/){$Rareness="Common"; $c_phase3_common++;}
	    $seenINdb++;
	  }

	  ##### 
	  if($Rareness eq "NA" && $seenINdb > 0){$Rareness="Rare"}
          # if there is any submission saying about sign
  
          if($clinVar > 0 ) { $clinvarness = "ClinVar" ; $c_clinsig_rare++;} 

	  if($Rareness =~ /Rare/){$c_SNV_match_1KG_AF_1++; $c_SNV_match_1KG++;}
	  elsif($Rareness =~ /Common/){$c_SNV_match_1KG++; $c_SNV_match_Common++;}
	  elsif($Rareness =~ /NA/){$c_SNV_nt_match++}
          
          # clinvar remove
          my $maf_upperThr = "hit";
          if($Phase3_af > $VAR{'clinVar_Max_MAF'} || $exac_AF > $VAR{'clinVar_Max_MAF'}) {
            $maf_upperThr = "notHit";
          }

	  if($Rareness eq "Rare" || $Rareness eq "NA" || ($clinvarness eq "ClinVar" && $maf_upperThr eq "hit" ) || ($VAR{'trioORsingle'} == 2 && $VAR{'D_model'}=~/somatic/i))
	  {
            push(@FilteringParameters, $Rareness); 
            push(@FilteringParameters, $clinvarness);

	    Filter_single_GT(@FilteringParameters), if($VAR{'trioORsingle'} == 0);
	    Filter_single_PatientOnly(@FilteringParameters), if($VAR{'trioORsingle'} == 1);
	    Filter_single_CancerSample(@FilteringParameters), if($VAR{'trioORsingle'} == 2);
	  }  
        }
      } 
    }
  }
}

close $IN;

print LOG "HEAD\tPatient_ID\tDisease_Model\tVariation\tAllRead_GT_MinDP\t";
print LOG "QUAL_GT_20\tPresentInDatabase\t";
print LOG "CommonInGnomAD_Exomes\tCommonInGnomAD_Genomes\tCommonInLC2\tCommonInLC3\tCommonInLC4\tClinVAR_Significant\t";
print LOG "CommonInPhase1\tCommonInPhase3\tCommonVariants\t";
print LOG "RareVariants\tPrivateVariants\t";
print LOG "ForFiltering\tFuntionalVariants\t";
print LOG "Denovo\tHomozygous\tHemizygous\tHeterozygous\n";

print LOG "SUM\t$VAR{'patientPIDs'}\t$VAR{'D_model'}\t$st_variation\t$c_Allread_gt10\t";
print LOG "$c_SNV_more_confidence_th\t$c_SNV_match_1KG\t";
print LOG "$c_ExAC_common\t$c_gnomAD_common\t$c_LC2_common\t$c_LC3_common\t$c_LC4_common\t$c_clinsig_rare\t";
print LOG "$c_phase1_common\t$c_phase3_common\t$c_SNV_match_Common\t";
print LOG "$c_SNV_match_1KG_AF_1\t$c_SNV_nt_match\t";
print LOG "$c_SNV_for_filtering\t$c_SNV_nsy_spl\t";

my @potential_D_models=qw/Denovo Homozygous Hemizygous Heterozygous/;
map{$c_SNV_filtered{$_} ? print LOG "$c_SNV_filtered{$_}\t" : print LOG "0\t"} @potential_D_models;
print LOG "\n";

close LOG;
close SUMMARY;

sub Filter_single_CancerSample
{
  $c_SNV_for_filtering++;
  
  my $D_model          = $_[0];
  my @query_FN         = @{$_[1]};
  my @pattern_FN       = @{$_[2]};
  my $variation        = $_[3];
  my $patientPIDs      = $_[4];
  my @DP_field         = @{$_[5]};
  my $st_D_model       = $_[6];
  my $st_variation     = $_[7];
  my @query_alt_AF     = @{$_[8]};
  my @BaseQ            = @{$_[9]};
  my $dbsnp            = $_[10];    
  my @query_dbSNP_1KG  = @{$_[11]};
  my @splitInputLine   = @{$_[12]};
  my $GQ_info          = $_[13];
  my $DP_AllAbove10or20 = $_[14];
  my $GQ_Format        = $_[15];
  my @FowreadRevrread  = @{$_[16]};
  my $rareness         = $_[17];
  my $clinvarness      = $_[18];

  my $varTag = join("_", @splitInputLine[0..1], @splitInputLine[3..4]);
  #print "CHECK\t$varTag\t$query_FN[0]=~/$pattern_FN[0]/ || $query_FN[1]=~/$pattern_FN[1]/\n";  
  if($query_FN[0]=~/$pattern_FN[0]/ || $query_FN[1]=~/$pattern_FN[1]/)
  { 
    $c_SNV_nsy_spl++;
    # Extracting Genotype information 
    # This will work for both Indel and SNVs
    my ($GQ, @GTs) = Check_GQ(\@DP_field, 2, $variation, \@query_alt_AF, $GQ_info, $splitInputLine[6], $GQ_Format);
    #print "CHECK\t$GTs[0]\t$GTs[1]\n";
    my @toPrintParameters = ($patientPIDs, $st_D_model, $st_variation, \@DP_field, \@query_FN, \@BaseQ, $dbsnp, \@query_dbSNP_1KG, \@query_alt_AF, \@splitInputLine, $DP_AllAbove10or20, \@GTs, \@FowreadRevrread, $rareness, $clinvarness);

    #if((($st_variation eq "SNVs" && $FowreadRevrread[0] eq "Yes" && $FowreadRevrread[1] eq "Yes") || $st_variation eq "Indel") || $D_model=~/Somatic/i)
    #{ 
      ###### Germline variants
      if($splitInputLine[6]!~/alleleBias/) {
        if($GTs[0] =~/HOMO_A|HOMO_A_lowLH/ && $GTs[1] =~/HOMO_A|HOMO_A_lowLH/)
        {
          #$c_SNV_filtered++;
          $toPrintParameters[1]='Homozygous';
          if($GTs[0]=~/lowLH/ || $GTs[1]=~/lowLH/){$toPrintParameters[1]='Homozygous_lowLH'}
          Print_Summary(@toPrintParameters);
        }
        if($GTs[0] =~/HETER|HETER_lowLH/ && $GTs[1] =~/HETER|HETER_lowLH/) # Heterozygous
        {         
          #$c_SNV_filtered++;
          $toPrintParameters[1]='Heterozygous';
          if($GTs[0]=~/lowLH/ || $GTs[1]=~/lowLH/){$toPrintParameters[1]='Heterozygous_lowLH'}
          Print_Summary(@toPrintParameters);
        }
        if($GTs[0]=~/HETER|HETER_lowLH/ && $GTs[1]=~/HOMO_A|HOMO_A_lowLH/)
        {			
          #$c_SNV_filtered++;
          $toPrintParameters[1]='LOR';
          if($GTs[0]=~/lowLH/ || $GTs[1]=~/lowLH/){$toPrintParameters[1]='LOR_lowLH'}
          Print_Summary(@toPrintParameters);
        }
        if($GTs[0]=~/HETER|HETER_lowLH/ && $GTs[1]=~/HOMO_R|HOMO_R_lowLH/)
        {			
          #$c_SNV_filtered++;
          $toPrintParameters[1]='LOA';
          if($GTs[0]=~/lowLH/ || $GTs[1]=~/lowLH/){$toPrintParameters[1]='LOA_lowLH'}
          Print_Summary(@toPrintParameters);
        }
      }
      ########### somatic
      if($GTs[0]=~/HOMO_R/ && $GTs[1]=~/HETER/ && $st_D_model=~/Somatic/)
      {
        #$c_SNV_filtered++;
        $toPrintParameters[1]='Heterozygous';
        if($GTs[0]=~/lowLH/ || $GTs[1]=~/lowLH/){$toPrintParameters[1]='Heterozygous_lowLH'}
        Print_Summary(@toPrintParameters);
      }
      if($GTs[0]=~/HOMO_R/ && $GTs[1]=~/HOMO/ && $st_D_model=~/Somatic/)
      {
        #$c_SNV_filtered++;
        $toPrintParameters[1]='Homozygous';
        if($GTs[0]=~/lowLH/ || $GTs[1]=~/lowLH/){$toPrintParameters[1]='Homozygous_lowLH'}
        Print_Summary(@toPrintParameters);
      }
      
      ## Denovo variants - Can take in alleleBias variants as well
      if($GTs[0]=~/HOMO_R/ && $GTs[1]=~/HOMO|HETER/ && $st_D_model=~/Germline/)
      {
        #$c_SNV_filtered++;
        $toPrintParameters[1]='Denovo';
        if($GTs[0]=~/lowLH/ || $GTs[1]=~/lowLH/){$toPrintParameters[1]='Denovo_lowLH'}
        Print_Summary(@toPrintParameters);
      }

    #}  
  }  
}


sub Filter_single_PatientOnly
{
  $c_SNV_for_filtering++;  
  my $D_model          = $_[0];
  my @query_FN         = @{$_[1]};
  my @pattern_FN       = @{$_[2]};
  my $variation        = $_[3];
  my $patientPIDs      = $_[4];
  my @DP_field         = @{$_[5]};
  my $st_D_model       = $_[6];
  my $st_variation     = $_[7];
  my @query_alt_AF     = @{$_[8]};
  my @BaseQ            = @{$_[9]};
  my $dbsnp            = $_[10];
  my @query_dbSNP_1KG  = @{$_[11]};
  my @splitInputLine   = @{$_[12]};
  my $GQ_info          = $_[13];
  my $DP_AllAbove10or20 = $_[14];
  my $GQ_Format        = $_[15];
  my @FowreadRevrread  = @{$_[16]};
  my $rareness         = $_[17];
  my $clinvarness      = $_[18];

  my $varTag = join("_", @splitInputLine[0..1], @splitInputLine[3..4]); 

  if($query_FN[0]=~/$pattern_FN[0]/ || $query_FN[1]=~/$pattern_FN[1]/)
  {
    $c_SNV_nsy_spl++;    
    # Extracting Genotype information 
    my ($childGT, $PL_HR, $PL_ALT, $PL_HA, $GQ, @GTs);
    my $DP;
    ($GQ, @GTs) = Check_GQ(\@DP_field, 1, $variation, \@query_alt_AF, $GQ_info, $splitInputLine[6], $GQ_Format);
    push(@GTs, ('NA', 'NA'));
    
    
    my @toPrintParameters = ($patientPIDs, $st_D_model, $st_variation, \@DP_field, \@query_FN, \@BaseQ, $dbsnp, \@query_dbSNP_1KG, \@query_alt_AF, \@splitInputLine, $DP_AllAbove10or20, \@GTs, \@FowreadRevrread, $rareness, $clinvarness);
    
    #if($GQ > $VAR{'genotypeQuality'} && (($st_variation eq "SNVs" && $FowreadRevrread[0] eq "Yes") || $st_variation eq "Indel"))
    if($GQ > $VAR{'genotypeQuality'})
    {
      if($GTs[0] =~/HOMO_A/)        
      {
        #$c_SNV_filtered++;          
        $toPrintParameters[1]="Homozygous";
        Print_Summary(@toPrintParameters);         
      }          
      # only Patient - Heterozygous
      # This will work for both Indel and SNVs
      if($GTs[0]=~/HETE|LOW_AF/)        
      {
        #$c_SNV_filtered++;         
        $toPrintParameters[1]="Heterozygous";
        Print_Summary(@toPrintParameters);          
      }            
      #count filtered models
      $c_SNV_filtered{$toPrintParameters[1]}++;
    } 
  }  
}

sub Filter_single_GT
{
  $c_SNV_for_filtering++;
    
  my $D_model          = $_[0];
  my @query_FN         = @{$_[1]};
  my @pattern_FN       = @{$_[2]};
  my $variation        = $_[3];
  my $patientPIDs      = $_[4];
  my @DP_field         = @{$_[5]};
  my $st_D_model       = $_[6];
  my $st_variation     = $_[7];
  my @query_alt_AF     = @{$_[8]};
  my @BaseQ            = @{$_[9]};
  my $dbsnp            = $_[10];    
  my @query_dbSNP_1KG  = @{$_[11]};
  my @splitInputLine   = @{$_[12]};
  my $GQ_info          = $_[13];
  my $DP_AllAbove10or20 = $_[14];
  my $GQ_Format        = $_[15];
  my @FowreadRevrread  = @{$_[16]};
  my $rareness         = $_[17];
  my $clinvarness      = $_[18];
 
  if($query_FN[0]=~/$pattern_FN[0]/ || $query_FN[1]=~/$pattern_FN[1]/)
  { 
    $c_SNV_nsy_spl++;
    
    # Extracting Genotype information 
    # This will work for both Indel and SNVs
    
    my ($GQ, @GTs) = Check_GQ(\@DP_field, 3, $variation, \@query_alt_AF, $GQ_info, $splitInputLine[6], $GQ_Format);
    my @toPrintParameters = ($patientPIDs, $st_D_model, $st_variation, \@DP_field, \@query_FN, \@BaseQ, $dbsnp, \@query_dbSNP_1KG, \@query_alt_AF, \@splitInputLine, $DP_AllAbove10or20, \@GTs, \@FowreadRevrread, $rareness, $clinvarness);
  
    $toPrintParameters[1]=""; 
    if($GQ > $VAR{'genotypeQuality'})
    {
      # Denovo 
      if($GTs[0]=~/HETER|HOMO_A|LOW_AF/ && $GTs[1] =~ /HOMO_R/ && $GTs[2] =~ /HOMO_R/)
      {          
        #$c_SNV_filtered++;
        $toPrintParameters[1]="Denovo";
        
        Print_Summary(@toPrintParameters);
      }
      # Homozygous
      if($GTs[0] =~/HOMO_A/ && $GTs[1] =~ /HETER|LOW_AF/ && $GTs[2]=~/HETER|LOW_AF/)
      {
        #$c_SNV_filtered++;
        $toPrintParameters[1]="Homozygous";
        Print_Summary(@toPrintParameters);
      }
      # Heterozygous
      if($GTs[0]=~/HETER|LOW_AF/ && (($GTs[1]=~/HETER|LOW_AF/ || $GTs[2]=~/HETER|LOW_AF/) || ($GTs[1]=~/HOMO_A/ || $GTs[2]=~/HOMO_A/)))
      {
        #$c_SNV_filtered++;
        $toPrintParameters[1]="Heterozygous";
        Print_Summary(@toPrintParameters);
      }
      # Hemizygous # Homo alternative in patient; Hetero in anyone of the parent other should be homo to reference ! 
      if($GTs[0]=~/HOMO_A/&& (($GTs[1]=~/HETER|LOW_AF/ && $GTs[2]=~/HOMO_R/) || ($GTs[1]=~/HOMO_R/ && $GTs[2]=~/HETER|LOW_AF/)))
      {
        #$c_SNV_filtered++;
        $toPrintParameters[1]="Hemizygous";
        Print_Summary(@toPrintParameters);
      }
      #count
      $c_SNV_filtered{$toPrintParameters[1]}++;
    }  
  }  
}

# To check whether the filtered variants are present in control parents
# For denovo - same variants with heterozygous in control parents
# For Homo - same variants with homozygous altervative in control parents 
# For Hetero - Not implemented
sub ParentControl
{
  my $st_variation   = $_[0];
  my $st_D_model     = $_[1];  
  my @splitInputLine = @{$_[2]};
  my $blackORlocal   = $_[3];
  
  my $inControl = "NoIn".$blackORlocal;
  
  # Control variant file created from parents variant calling from RaKi
  my ($allParents);
  my $key = join("_", @splitInputLine[0..1], @splitInputLine[3..4]);
  
  my $childorParent;  
  if($blackORlocal    =~/Control_1/){$allParents = $P_RAKI_UGD{$key}}
  elsif($blackORlocal =~/Control_2/){$allParents = $P_HIPO_EURO{$key}}
  elsif($blackORlocal =~/BlackList/){$allParents = $BlackList{$key}} 

  my @allParentDPinfo =split(/;/, $allParents), if($allParents!~/^$/);
  
  my (@allVAF, @allFilter, $ControlDPtoPrint);
  
  if($st_variation=~/SNVs/)
  {
    @allVAF = map {CalculateVAF($_)} @allParentDPinfo;    
    $ControlDPtoPrint = $allParents;
    
    if($st_D_model=~/Denovo|Hetero|Somatic/i)
    {
      @allFilter= map {($_ >= $VAR{'MinVAF_HT'} && $_ <= $VAR{'MaxVAF_HT'}) ? "YES" : "NO"} @allVAF;      
    }
    if($st_D_model=~/Homo|Hemi|LOH/i)
    {
      @allFilter= map {($_> $VAR{'MaxVAF_HoAlt'}) ? "YES" : "NO"} @allVAF;      
    }
  }  
  if($st_variation =~/Indel/)
  {
    $ControlDPtoPrint = $allParents;
    
    if($st_D_model=~/Denovo|Hetero|Somatic/i)
    {
      @allFilter=map{(/^0\/1|^1\/0/) ? "YES": "NO"} @allParentDPinfo;
    }
    if($st_D_model=~/Homo|Hemi|LOH/i)
    {
      @allFilter=map{(/^1\/1/) ? "YES": "NO"} @allParentDPinfo;      
    }
  }
  my $totalCount = 0;
  my $yesCount = 0;
  
  $totalCount = $#allFilter+1;
    
  map{/YES/? $yesCount++:"NA"} @allFilter;
    
  if(grep(/YES/, @allFilter))
  #if($totalCount >= 1)
  {
    $inControl = "YesIn".$blackORlocal;
  }
  # returns the binary results, all dp fields, total number and yes numbers 
  return($inControl, $ControlDPtoPrint, $totalCount, $yesCount);
}

# Calulated the Variant Allele frequency from the INFO DP field
sub CalculateVAF
{
  my $DP_field = $_[0];  
   
  my ($fr, $rr, $fnr, $rnr);
  
  ($fr, $rr, $fnr, $rnr)=$DP_field=~/DP\d=(\d+),(\d+),(\d+),(\d+)/; 
  
  my $DP = $fr+$rr+$fnr+$rnr;
  my $alt_AF = -1;
    
  if($DP > $VAR{'MinReads'})
  {
    $alt_AF = ($fnr+$rnr)/($DP);
  }
  
  return($alt_AF);
}

# Unique array
sub uniq 
{
  return keys %{{ map { $_ => 1 } @_ }};
}

# Check for the presence of variants from BlackList file provied by Matthias
# Similar to control parent check
# For Heterozygous - Not implemented
sub CheckBlackList
{  
  my $st_variation   = $_[0];
  my $st_disease     = $_[1];
  my $patientPID     = $_[2];
  my @splitInputLine = @{$_[3]};
  
  my $analysisPath = $VAR{'CLUSTER_EO'};
  my $tempFile = "$analysisPath/temp_$patientPID.$st_disease.$st_variation.vcf";
  
  open(TMP, ">$tempFile") || die "Can't create the temp file $tempFile\n";
  
  my @newInputLines = @splitInputLine[0..7];
  
  # For now lets use this (its FAKE header)
  my $VCFheader = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO";
  print TMP $VCFheader, "\n";
  print TMP join("\t", @newInputLines), "\n";
  close TMP;
  
  my $BINannotate_vcf = $VAR{'SCRIPT_ANNOTATE_VCF'};
  my $blackListVcf = $VAR{'BLACK_LIST_VCF'};
  chomp(my $blackList=`perl $BINannotate_vcf -b $blackListVcf --bFileType vcf -a $tempFile --aFileType vcf --columnname BlackList --bAdditionalColumns 5|grep -v '#'`);
 
  my @split_blackList=split(/\t/, $blackList);
  # Removing the tmp file;
  `rm $tempFile`;
  
  my $colInfo = $split_blackList[$#split_blackList];
  my (@ReadsVAF, $DP4Reads);
  if($colInfo=~/^MATCH/)
  {  
    ($DP4Reads) = $colInfo=~/DP4READS=(.*)$/;  
    my @splitDP4Reads = split(/;/, $DP4Reads);
    @ReadsVAF = map {CalculateVAF($_)} @splitDP4Reads;
  }
  
  my @Filter = ();
  
  if($st_disease =~/Denovo|Hetero|Somatic/)
  {    
    @Filter= map {($_ >= $VAR{'MinVAF_HT'} && $_ <= $VAR{'MaxVAF_HT'}) ? "YES" : "NO"} @ReadsVAF;    
  }
  if($st_disease =~/Homo|Hemi|LOH/)
  {   
    @Filter= map {($_ > $VAR{'MaxVAF_HoAlt'}) ? "YES" : "NO"} @ReadsVAF;    
  }
  
  my $inBlackList = "NoInBlackList";
  my $totalCount = 0;
  my $yesCount = 0;
  $totalCount = $#Filter+1;
  map{/YES/ ? $yesCount++ : "NA"} @Filter;

  if(grep(/YES/, @Filter))
  #if($totalCount >=1 )
  {
    $inBlackList = "YesInBlackList";    
  }
  return($DP4Reads, $inBlackList, $totalCount, $yesCount);
}

# Check for Genotype quality and Genotypes from DP fields
sub Check_GQ
{
  my @DP_field          = @{$_[0]};  
  my $number_of_samples = $_[1];
  my $variation         = $_[2];
  my @query_alt_AF      = @{$_[3]};
  my $GQ_INFO           = $_[4];
  my $IndelFilter       = $_[5];
  my $GQ_Format         = $_[6];

  my $PL_HARDFilter     = -4.60517;
  #my $PL_HARDFilter     = 0;
  my $GT_quality        = $VAR{'genotypeQuality'};
  
  my $Filter_PASS;
  
  if($VAR{'trioORsingle'} == 2 && $VAR{'D_model'}=~/somatic/i)
  {
    $Filter_PASS       = ".";  
  }
  elsif($VAR{'trioORsingle'} == 2) 
  {
    $Filter_PASS       = "PASS|alleleBias|QD;alleleBias";
  }
  else
  {
    $Filter_PASS       = "PASS";
  }

  my $GQ;
#  my ($iGT, $iGQ, $iPL, $iNV, $iDP);
  my $DP_HLHS;
  my ($PL_R, $PL_HE, $PL_A);
  my @GTs = ();
  
  #my @ssFormat = split(/:/, $GQ_Format);
  #my @ssInfo   = split(/:/, $GQ_INFO);

  #for(my $i=0; $i<=$#ssFormat; $i++)
  #{
  #  if($ssFormat[$i] eq "GT"){$iGT=$ssInfo[$i]}
  #  if($ssFormat[$i] eq "GQ"){$iGQ=$ssInfo[$i]}
  #  if($ssFormat[$i] eq "PL" || $ssFormat[$i] eq "GL"){$iPL=$ssInfo[$i]}    
  #}  

  my $count=0;
  if(($variation =~/Indel|SNVs/ && $VAR{'trioORsingle'} < 2)|| ($variation=~/SNVs|Indel/ && $VAR{'trioORsingle'} == 2 && $VAR{'D_model'} !~/Somatic/i ))
  #if($variation =~/Indel/ || ($variation=~/SNVs/ && $VAR{'D_model'}!~/Somatic/ && $VAR{'trioORsingle'} == 2))
  {
    for(my $k=0; $k<=$#DP_field; $k++)
    {
      my $GT = 'NA';
      my ($iGT, $iGQ, $iPL, $iNV, $iDP);
      if($DP_field[$k]!~/NA/)
      {        
        my @ssFormat = split(/:/, $GQ_Format);
        my @ssInfo   = split(/:/, $DP_field[$k]);
  
        for(my $i=0; $i<=$#ssFormat; $i++)
        {
          if($ssFormat[$i] eq "GT"){$iGT=$ssInfo[$i]}
          if($ssFormat[$i] eq "GQ"){$iGQ=$ssInfo[$i]}
          if($ssFormat[$i] eq "PL" || $ssFormat[$i] eq "GL"){$iPL=$ssInfo[$i]}
          if($ssFormat[$i] =~ /NV|AD/){$iNV=$ssInfo[$i]}
          if($ssFormat[$i] =~/NR|DP/){$iDP=$ssInfo[$i]}
        }
        if($k==0){$DP_HLHS = $iDP}
  
        #if($field eq 'NA'){next;}
        #($PL_R, $PL_HE, $PL_A, $GQ)=$field=~/:(.?\d+\.?\d+?),(.?\d+\.?\d+?),(.?\d+\.?\d+?):\d+:(\d+):\d+:\d+$/;
        ($PL_R, $PL_HE, $PL_A)=split(/,/, $iPL);
        #($PL_R, $PL_HE, $PL_A)=$iPL=~/(.?\d+\.?\d+?),(.?\d+\.?\d+?),(.?\d+\.?\d+?)/;      
        
        #my @split_field=split(/:/, $field);
        
  
        # Basic platypus filtering. Since we dont use the PASS in filter in platypus anymore. 
        # Tumor samples can go beyond PASS filter
        if($IndelFilter =~/$Filter_PASS/)
        {
          if($iGT=~/0[\/|\|]1|1[\/|\|]0/)
          {
           if($PL_R  < $PL_HARDFilter && $PL_A < $PL_HARDFilter){$GT="HETER";}else{$GT="HETER_lowLH";} 
  #         $GT="HETER";
          }
          elsif($iGT=~/0[\/|\|]0/)
          {
            if($PL_HE < $PL_HARDFilter && $PL_A < $PL_HARDFilter){$GT="HOMO_R";}else{$GT="HOMO_R_lowLH"}
  #           $GT="HOMO_R";
          } 
          elsif($iGT=~/1[\/|\|]1|^1$/)
          {
            if($PL_HE < $PL_HARDFilter && $PL_R < $PL_HARDFilter){$GT="HOMO_A";}else{$GT="HOMO_A_lowLH"}
  #          $GT="HOMO_A";
          }
        }
      }  
      push(@GTs, $GT);
      $count++, if($iGQ > $GT_quality);        
    }
    # Filter passes if all 3 samples (TRIO samples) or 1 sample (only patient) have better genotype quality    
   if($count==$number_of_samples){$GQ=$GT_quality+1}; 
  }
  elsif($variation=~/SNVs|Indel/ && $VAR{'trioORsingle'}  == 2 && $VAR{'D_model'} =~/Somatic/i)
  #elsif($variation =~/SNVs/)
  {
    #($GQ)=$GQ_INFO=~/:\d+:\d+:(\d+)$/;
    ($GQ)=$GQ_INFO=~/:(\d+)$/;
    #$GQ = $iGQ;
    foreach my $AF(@query_alt_AF)   
    {
      my $GT="NA";
      if($AF!~/NA/)
      {
        if($AF >= $VAR{'MinVAF_HT'} && $AF <= $VAR{'MaxVAF_HT'}){ $GT="HETER"}
        elsif($AF > 0 && $AF < $VAR{'MaxVAF_HoRef'}){ $GT="LOW_AF"}
        elsif($AF == 0){$GT="HOMO_R"}
        elsif($AF > $VAR{'MaxVAF_HoAlt'}){ $GT="HOMO_A"}
      }   
      push(@GTs, $GT);
    }
  }
  return($GQ, @GTs);
}

# Printing the summary 
# Which also calls for checkBlackList, GeneticAssocation, ParentControl
# Calling the above functions here makes more sense
sub Print_Summary
{
  my $patientPIDs     = $_[0];
  my $st_D_model      = $_[1];
  my $st_variation    = $_[2];  
  my @DP_field        = @{$_[3]};
  my @query_FN        = @{$_[4]};
  my @BaseQ           = @{$_[5]};
  my $dbsnp           = $_[6];  
  my @query_dbSNP_1KG = @{$_[7]};
  my @query_alt_AF    = @{$_[8]};
  my @splitInputLine  = @{$_[9]};
  my $DP_AllAbove10or20 = $_[10];
  my @GTs               = @{$_[11]};
  my @FowreadRevrread   = @{$_[12]};
  my $rareness          = $_[13];
  my $clinvarness       = $_[14];

  
  my $child_DP_field  = $DP_field[0];
  my $mother_DP_field = $DP_field[1];
  my $father_DP_field = $DP_field[2];

  my ($child_DP, $child_alt_AF, $child_BaseQ)    = Extract_DP($child_DP_field, $VAR{'variation'}, $VAR{'D_model'}, $VAR{'trioORsingle'});

  my $child_BaseQ  = $BaseQ[0];
  my $mother_BaseQ = $BaseQ[1];
  my $father_BaseQ = $BaseQ[2];
   
  
  my ($blackListDPtoPrint, $geneticAssociation, $inBlackList);
  my ($BLTotalCount, $BLYesCount);
  my ($inControl_RU, $controlDPtoPrint_RU, $PCTotalCount_RU, $PCYesCount_RU);
  my ($inControl_HE, $controlDPtoPrint_HE, $PCTotalCount_HE, $PCYesCount_HE);
  my ($IndividualLikelihoodInControl, $mutationRateInGene);
  
  my $varTag = join("_", @splitInputLine[0..1], @splitInputLine[3..4]);
  
  my $geneName = $splitInputLine[$gene_VcfCol];
  my $transcriptData = $splitInputLine[$Transcript_VcfCol];
  #my $CADD_Phred = GetCADDScore($varTag);
  #my $promoterGene =InPromoter($varTag);
  
  # This pattern will get only the first gene, in case of intergenic variants
   #$geneName=~s/\(dist=(\d+|NONE)\)//g;
   $geneName=~s/\(.*?\)//g;
#  if($geneName=~/(\w+-?\.?\w+)\(dist=\d+|NONE\),(\w+-?\.\w+)\(dist=\d+|NONE\)/){$geneName="$1".","."$2"}
#  elsif($geneName=~/(\w+-?\.?\w+);(\w+-?\.?\w+)/){$geneName="$1".","."$2"}
#  elsif($geneName=~/(\w+-?\.?\w+)[;|,|\(]?/){$geneName=$1}
  
  ### Checking with the local controls##################################
  if($VAR{'RUN_BLACK_LIST'} == 1)
  {
    #($blackList, $inBlackList, $BLTotalCount, $BLYesCount) = CheckBlackList($st_variation, $st_D_model, $patientPIDs, \@splitInputLine);
    ($inBlackList, $blackListDPtoPrint, $BLTotalCount, $BLYesCount) = ParentControl($st_variation, $st_D_model, \@splitInputLine, 'BlackList');
  }  
  
  # Local exome controls from UniExome project
  if($VAR{'RUN_LOCAL_CONTROL'} == 1)
  {
    ($inControl_RU, $controlDPtoPrint_RU, $PCTotalCount_RU, $PCYesCount_RU) = ParentControl($st_variation, $st_D_model, \@splitInputLine, 'Control_1');
    ($inControl_HE, $controlDPtoPrint_HE, $PCTotalCount_HE, $PCYesCount_HE) = ParentControl($st_variation, $st_D_model, \@splitInputLine, 'Control_2');
  } 
   
  my $inBC = "NA";
  if($inBlackList=~/NoInBlackList/ && $inControl_HE=~/NoInControl/ && $inControl_RU=~/NoInControl/)
  {
    $c_Finalfiltered++;
    $inBC = "NoInBC";
  }
  elsif($inBlackList=~/YesInBlackList/ || $inControl_HE=~/YesInControl/ || $inControl_RU=~/YesInControl/)
  {
   $inBC = "YesInBC"; 
  }
  
  # NA when CONTROL_INFO is empty
  if($controlDPtoPrint_RU=~/^$/){$controlDPtoPrint_RU = "NA"}
  if($controlDPtoPrint_HE =~/^$/){$controlDPtoPrint_HE = "NA"}
  if($blackListDPtoPrint  =~/^$/){$blackListDPtoPrint  = "NA"}
  
  #Start printing 
  my $key = join("_", $splitInputLine[0],$splitInputLine[1],$splitInputLine[3],$splitInputLine[4]);

  if($VAR{'toPrintSummaryORVCF'} == 0 )
  {
    print SUMMARY "var_$CaseControl{$patientPIDs}\t$patientGroup{$patientPIDs}\t$patientPIDs\t$VAR{'D_model'}\t$st_D_model\t$st_variation\t$dbsnp\t";
    print SUMMARY "$splitInputLine[0]\t$splitInputLine[1]\t$splitInputLine[3]\t$splitInputLine[4]\t";
    print SUMMARY "$geneName\t$child_DP\t$child_DP_field\t$mother_DP_field\t$father_DP_field\t";
    
    #print SUMMARY "$child_BaseQ\t$mother_BaseQ\t$father_BaseQ\t";
    print SUMMARY join ("\t", @query_alt_AF), "\t";
    print SUMMARY "$GTs[0]\t$GTs[1]\t$GTs[2]\t";

    print SUMMARY "$query_FN[0]\t$query_FN[1]\t";
    print SUMMARY "$transcriptData\t$query_dbSNP_1KG[0]\t$query_dbSNP_1KG[1]\t";

    # ExAC - gnomAD exomes 
    if(defined $ExACdata{$varTag}{'MAF'}){print SUMMARY "$ExACdata{$varTag}{'MAF'}\t$ExACdata{$varTag}{'homo'}\t$ExACdata{$varTag}{'info'}\t";}
    else{print SUMMARY ".\t.\t.\t";}

    #gnomAD genomes
    if(defined $gnomADdata{$varTag}{'MAF'}){print SUMMARY "$gnomADdata{$varTag}{'MAF'}\t$gnomADdata{$varTag}{'homo'}\t$gnomADdata{$varTag}{'info'}\t";}
    else{print SUMMARY ".\t.\t.\t";}

    # local control 2 filtering
    if(defined $LocalControl_2{$varTag}{'AF'}){print SUMMARY "$LocalControl_2{$varTag}{'AF'}\t"}
    else {print SUMMARY ".\t"}

    # local control 3 filtering
    if(defined $LocalControl_3{$varTag}){print SUMMARY "$LocalControl_3{$varTag}\t"}
    else {print SUMMARY ".\t"}

    # Phase 3 filtering
    if(defined $Phase3{$varTag}){print SUMMARY "$Phase3{$varTag}\t$splitInputLine[$phase3_VcfCol]\t"; }
    else { print SUMMARY ".\t.\t";}

    # local control 4 filtering
    if(defined $LocalControl_4{$varTag}){print SUMMARY "$LocalControl_4{$varTag}\t"}
    else{print SUMMARY ".\t";}

    print SUMMARY "$rareness\t$clinvarness\t";  
    #print SUMMARY "$BLYesCount/$BLTotalCount\t$inBlackList\t";
    #print SUMMARY "$PCYesCount_RU/$PCTotalCount_RU\t$inControl_RU\t";
    #print SUMMARY "$PCYesCount_HE/$PCTotalCount_HE\t$inControl_HE\t";
    #print SUMMARY "$inBC\t";
    
    if(defined $GADInfo{$geneName})
    {
      my @uniqGADInfo = uniq(@{$GADInfo{$geneName}});
      print SUMMARY join(", ", @uniqGADInfo), "\t";
    }
    else{print SUMMARY ".\t"}
    
    if(defined $DisGeNETInfo{$geneName}){print SUMMARY join(", ", @{$DisGeNETInfo{$geneName}}), "\t"}else{print SUMMARY ".\t"}
    if(defined $OMIM_Number{$geneName}){ print SUMMARY "$OMIM_Number{$geneName}\t"}else{print SUMMARY ".\t"}
    if(defined $OMIM_Disorder{$geneName}){ print SUMMARY "$OMIM_Disorder{$geneName}\t"}else{print SUMMARY ".\t"}
    if(defined $clinVarInfo{$varTag}){print SUMMARY "$clinVarInfo{$varTag}\t"}else{print SUMMARY ".\t.\t"}

    if(defined $RVIS_LC{$geneName}){print SUMMARY "$LocalControlMutation{$geneName}\t$RVIS_LC{$geneName}\t"}
    else{print SUMMARY ".\t.\t.\t"};
    if(defined $RVIS_PG{$geneName}){print SUMMARY "$RVIS_PG{$geneName}\t"}else{print SUMMARY ".\t"};
    if(defined $RVIS_ExAC{$geneName}){print SUMMARY "$RVIS_ExAC{$geneName}\t"}else{print SUMMARY ".\t"};

    if(defined $GeneDesc{$geneName}){print SUMMARY "$GeneDesc{$geneName}\t"}else{print SUMMARY ".\t"}
    if(defined $GOinfo{$geneName}){print SUMMARY join (", ", @{$GOinfo{$geneName}}), "\t"}else{print SUMMARY ".\t"}
    if(defined $bioGPS{$geneName}){print SUMMARY "$bioGPS{$geneName}\t"}else{print SUMMARY ".\t"}
    if(defined $SCLinfo{$geneName}){print SUMMARY "$SCLinfo{$geneName}\t"}else{print SUMMARY ".\t"}
    if(defined $typeInfo{$geneName}){print SUMMARY "$typeInfo{$geneName}"}else{print SUMMARY "."}
    #print SUMMARY "\t$CADD_Phred\t$promoterGene\n"; 
    #print SUMMARY "\t", join("\t", @FowreadRevrread);
    print SUMMARY "\n";
  }
  elsif($VAR{'toPrintSummaryORVCF'} == 1)  
  {
    print SUMMARY join("\t", @splitInputLine), "\n";
  }
}

# Warning and Die message
sub USAGE
{
    print <<USAGE_MESSAGE;
Arguments are not proper. Please check for correct usage.

USAGE: perl Filtering_VCF_file.pl Arguments
    1) Config file,with stable variants
    2) Input VCF file 
    3) Summary file/Output file
    4) Patient PIDS 
    5) Variation (snp==0 or indel==1) 
    6) Disease model (denovo/Homo-/Hetero-/Hemi- zygous)     
    7) Indel calling (Platypus==2)    
USAGE_MESSAGE
    die;
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


# Notes on Extract_DP
# Arguments:
#   1. DP field, DP4 or DP5
#   2. Dindel or Platypus
#   3. SNVs or Indels

sub Extract_DP
{
  my $DP_field = $_[0];  
  my $variation = $_[1];
  my $D_model = $_[2];
  my $trioORsingle =$_[3];
  
  my $DPall = 0;
  my $DP = 0;
  my $alt_AF = 0;
  my $BaseQ = 0;
  my $FowreadRevrread=".";
  
  # Mpileup SNV
  if($variation=~/SNVs|Indel/ && $trioORsingle  == 2 && $D_model =~/Somatic/i) 
  {
    ($DPall)=$DP_field=~/DP=(\d+);/;
    
    my ($fr, $rr, $fnr, $rnr);
    ($fr, $rr, $fnr, $rnr)=$DP_field=~/DP\d=(\d+),(\d+),(\d+),(\d+)/;    

    $DP = $fr+$rr+$fnr+$rnr;
  
    $alt_AF = ($fnr+$rnr)/($DP), if($DP>0);
  
    $BaseQ = sprintf("%.2f", $DP/$DPall), if($DPall>0);    
  
    if($fnr > 0 && $rnr > 0) {$FowreadRevrread="Yes"}else{$FowreadRevrread="No"}
  }
  # Platypus Indel and SNVs
  if(($variation =~/Indel|SNVs/ && $trioORsingle < 2)|| ($variation=~/SNVs|Indel/ && $trioORsingle == 2 && $D_model !~/Somatic/i ))
  {
    my $NonRef;
    # Multiple alternative allele
    if($DP_field=~/,/)
    {
      my @split_DP_field = split(/:/, $DP_field);
      
      my @splitNR = split(/,/, $split_DP_field[4]);
      my @splitNV = split(/,/, $split_DP_field[5]);
      
      my ($sumNR, $sumNV) = (0,0);
      $sumNR += $_ foreach @splitNR;
      $sumNV += $_ foreach @splitNV;
      
      if($#splitNR >=0){$DP = $sumNR/($#splitNR+1)}
      if($#splitNR >=0){$NonRef = $sumNV/($#splitNR+1)}
      $alt_AF=$NonRef/$DP, if($DP>0);
    }
    # Single alternative allele
    else
    {
      ($DP, $NonRef)=$DP_field=~/:(\d+):(\d+)$/;
      $alt_AF = $NonRef/$DP, if($DP>0);
    }    
  }
  return($DP, $alt_AF, $BaseQ, $FowreadRevrread);
}

sub GetCADDScore
{
  my @var=split("_", $_[0]);
  my $databaseSNVs = $VAR{'CADD_SNVs'};
  my $databaseIndel = $VAR{'CADD_Indel'};
  
  my @cadds="";
  #my @cadds = `tabix $databaseSNVs $var[0]:$var[1]-$var[1]; tabix $databaseIndel $var[0]:$var[1]-$var[1]`;
  my $CADD_phred="NoCADD";

  foreach my $cadd(@cadds)
  {
    chomp $cadd;
    my @caddScore = split(/\t/, $cadd);
    if($caddScore[0]==$var[0] && $caddScore[1]==$var[1] && $caddScore[2] eq $var[2] && $caddScore[3] eq $var[3])
    {
      $CADD_phred=$caddScore[5];
    }
  }
  return($CADD_phred);
}

## Variant in promoter  ? 
sub InPromoter
{
  my @var=split("_", $_[0]);
  my $promoterBED=$VAR{'PromoterFile'};
  my @promoters="";
  my $genePromoter="NotInPromoter";

  if($promoters[0]!~/^$/)
  {
    $genePromoter=join(";", @promoters);
  }
  return($genePromoter);
}

## ExAC population frequency
sub ExAC_sub
{
  my $ExAC_info = $_[0];
  my $maf_thr = $_[1];
  my $all_pop_maf=0;
  my $AF= ".";
  my $homo = ".";

  my $all_pop_maf_tag;

  if($ExAC_info=~/MATCH=exact/)
  {
    if($VAR{'runPopulationBasedFilter'}==1) {
      ($AF, $homo)=$ExAC_info=~/$VAR{'subPopulation'}_AF=(\d\.\d+);$VAR{'subPopulation'}_HOMO=(\d+);/;
      my @ss=split(/;/, $ExAC_info);

      ### checks for each populations whether AF greater 
      foreach my $ele(@ss)
      {
        if($ele=~/_AF=/)
        {
          my ($af)=$ele=~/_AF=(\d\.\d+(e-\d+)?)/;
          if($af > $maf_thr){$all_pop_maf++}
        }
      }
    }
    else
    {
      ($AF)=$ExAC_info=~/;AF=(\d\.\d+(e-\d+)?)/;
      ($homo) = $ExAC_info=~/;HOMO=(\d+);/;
       if($AF > $maf_thr){$all_pop_maf_tag="Common"}else{$all_pop_maf_tag="Rare"}
    }

    #if($AF > $maf_thr){$all_pop_maf_tag="Common"}else{$all_pop_maf_tag="Rare"}
  }
  return($AF, $homo, $all_pop_maf_tag);
}

#---------------------------------------
sub LC_sub
{
  my $LC_info = $_[0];
  my $AF = ".";
  if($LC_info=~/&/) {
    my @split_lc_info = split(/&/, $LC_info);
    foreach my $lc_info (@split_lc_info) {
      if($lc_info=~/MATCH=exact/)
      {
        my $AF_temp=".";
        ($AF_temp)=$lc_info=~/AF=(\d\.\d+)/;
        $AF = $AF_temp, if($AF_temp > $AF);
      }
    }
  }
  else{
    ($AF)=$LC_info=~/AF=(\d\.\d+)/;
  }
  return($AF);
}

#---------------------------------------
sub Phase3_sub
{
  my $phase3_info = $_[0];
  my $maf_thr = $_[1];
  my $all_pop_maf=0;
  my $AF = ".";

  my %samplePop = ('AFR' => 'AFR',
                'AMR' => 'AMR',
                'EAS' => 'EAS',
                'NFE' => 'EUR',
                'ASA' => 'SAS');

  if($phase3_info=~/MATCH=exact/)
  {
    ($AF)=$phase3_info=~/$samplePop{$VAR{'subPopulation'}}_AF=(\d\.?\d{0,})/;
    my @ss=split(/;/, $phase3_info);
    foreach my $ele(@ss)
    {
      if($ele=~/_AF=/)
      {
        my ($af)=$ele=~/_AF=(\d\.\d+)/;
        if($af > $maf_thr){$all_pop_maf++}
      }
    }
  }
  my $all_pop_maf_tag;
  if($VAR{'runPopulationBasedFilter'}==0)
  {
    #if($all_pop_maf > 0){$all_pop_maf_tag="Common"}else{$all_pop_maf_tag="Rare"}   
    ($AF) = $phase3_info=~/;AF=(\d\.?\d{0,})/;
    if($AF > $maf_thr){$all_pop_maf_tag="Common"}else{$all_pop_maf_tag="Rare"}
  }
  else
  {
    if($AF > $maf_thr){$all_pop_maf_tag="Common"}else{$all_pop_maf_tag="Rare"}
  }
  return($AF, $all_pop_maf_tag);
}

## ClinVar CLNSIG tag, 
### Variant Clinical Significance, 0 - Uncertain significance, 1 - not provided, 2 - Benign, 3 - Likely benign, 4 - Likely pathogenic, 5 - Pathogenic, 6 - drug response, 7 - histocompatibility, 255 - other"

sub clinSig 
{
  my $clinvar_info = $_[0];
  my $significant  = 0;
  my @a_clnsig       = ();

  if($clinvar_info =~ /MATCH=exact/) {
    my @multipleMatched =split(/\&/, $clinvar_info);
    foreach my $match(@multipleMatched) {

      my ($clnsig) = $match =~ /CLNSIG=(.*);CLNDSDB=/;
      push(@a_clnsig, $clnsig);

      my @clnsig_array = split(/[\||,]/, $clnsig);
 
      foreach my $submission(@clnsig_array) {
        if($submission != 2 && $submission != 3 && $submission != 255) {
          $significant += 1; 
        }
        else {
          $significant -= 1;
        }
      }
    }
  }
  my $join_a_clnsig = join("&", @a_clnsig);
  if($join_a_clnsig =~/^$/) {$join_a_clnsig = "."}

  return($significant, $join_a_clnsig);
}
