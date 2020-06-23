#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;

#########################################################################################
######## add all entries from dbNSFP to a vcf file with (non-synonymous) SNV     ######## 
######## by Mehmet Ates and Barbara Hutter, September 2013                       ######## 
#########################################################################################

if (@ARGV < 1)
{
  die "$0 needs at least 1 argument: \n\t1) vcf file with SNVs to add dbNSFP entries to - optional \n\t2) path to the bgzipped, tabix indexed dbNSFP database file (default: ./CADD/whole_genome_SNVs.tsv.gz)\n";
}

my $databasePATH="./CADD/";
chomp(my $inputfile = $ARGV[0]);
chomp(my $database = "$databasePATH/$ARGV[1]");

my $cmd; # tabix call;
unless (-r $database)
{
  die "could not read the dbNSFP database $database!";
}
# get header of dbNSFP
my $header_cmd = "tabix -h $database 1:0-0 |";  # pipe
open(HEAD, $header_cmd);  # laesst sich wie ein file handle behandeln!
my $header = <HEAD>;  # hat nur eine Zeile
close HEAD;
my @databasespalten = split ("\t", $header);
my $cols = @databasespalten;  # 52 fields in dbNSFP2.0

# go over the SNV file
my $EING;

if($inputfile=~/.gz$/)
{
  open ($EING, 'gzip -dc '. $inputfile."| ") or die "Could not open SNV file $inputfile!\n";
}
else
{
  open ($EING, $inputfile) or die "Could not open SNV file $inputfile!\n";
}

my $kommentare = 1;  # vcf header lines before "#CHROM"
my @spalten = ();
my $alteheader = "";
my $counter = 0;
my $present = 0;
my $zeile = "";
my $chro = "";
my $ref = "";
my $alt = "";
my $datazeile = "";
while (<$EING>)
{
  $zeile = $_;
  @spalten = split ("\t", $_);
  $alteheader = $_;
  if ($spalten[0] ne "#CHROM" && $kommentare == 1)  # vcf header lines before "#CHROM"
  {
    print;
    next;
  }
  elsif ($kommentare == 1 && $spalten[0] eq "#CHROM")  # last vcf header
  {
    $kommentare = 0;
    chomp ($alteheader);
    print "$alteheader\t$header";  # add dbNSFP header columns to last vcf header
    next;
  }
  ############ get fitting entries by tabix  
  if ($spalten[0] =~ /([\dXY]+)/ )  # get plain chromosome name (without "chr" or suffix)
  {
    $counter++;
    chomp $zeile;
    $chro = $1;
    $ref = uc $spalten[3];
    $alt = uc $spalten[4];
    $cmd = "tabix $database $chro:$spalten[1]-$spalten[1] |";
    open (RES, $cmd);
    my $gefunden = 0;
    my $existiert = 0;
    while (<RES>)
    {
      $existiert = 1;
      $datazeile = $_;
      @databasespalten = split /\t/, $_;
      # compare ref and alt allele of SNV file with dbNSFP to get the fitting one
      # there might be > 1 entry, such as for BRAF V600E (7:140453136-140453136)
      # the second is X>R, seems to be a ENSEMBL transcript frameshift: codon TGA instead of GTG
      # this one only has a mutation assessor score, no other, also has no Uniprot_acc
      # or FGFR1 8:38272308, same thing
      # MAP4K4 (2:102460600) has no Uniprot_acc in both, but here the first is the weird one. P>L is wrong, it's a R>X
      # and such entries are not "nonsynonymous" and don't get scores except LRT_pred => this is the one to look at,
      # it never has an entry for the "weird" things, but has one for stopgain
      if($alt=~/,/)
      {
        my @split_alt=split(/,/, $alt);
        my $high_cadds_alt;
        my $lastHigh=0;

        foreach my $alts(@split_alt)
        {
          if($ref eq $databasespalten[2] && $alts eq $databasespalten[3])
          {
            if($databasespalten[5] > $lastHigh){$lastHigh=$databasespalten[5]; $high_cadds_alt=$alts}
          }
          
        }
        print "$zeile\t$databasespalten[0]\t$databasespalten[1]\t$databasespalten[2]\t$databasespalten[3]\t$high_cadds_alt\t$lastHigh\n";
        $gefunden = 1;
        $present++;
      }
      else
      {
        if ($ref eq $databasespalten[2] && $alt eq $databasespalten[3])
        {
          print "$zeile\t$datazeile";
          $gefunden = 1;
          $present++;
        }
      }
    }
    close RES;
    if ($existiert == 0 || $gefunden == 0)  # there is either no entry or no matching one => fill the fields with \t.
    {
      print $zeile, "\tNoCADD" x $cols, "\n";
    }
  }
}
close $EING;
print STDERR "From $counter entries in SNV file, $present have a CADD entry.\n";
exit;
