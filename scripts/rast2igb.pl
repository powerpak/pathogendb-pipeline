#!/usr/bin/perl

# 08.01.2015 06:07:03 PST
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GLOBALS
my $sSvrRetrieveJob = 'svr_retrieve_RAST_job';
my $sGffToBed       = 'gff2bed.pl';
my $sFaToTwoBit     = 'faToTwoBit';

# GET PARAMETERS
my $sHelp        = 0;
my $sRastUser    = '';
my $sRastPass    = '';
my $nRastJobID   = '';
my $sGenomeName  = '';
my $sIGBdir      = '';
my $res = GetOptions("help!"    => \$sHelp,
                     "user=s"   => \$sRastUser,
                     "pass=s"   => \$sRastPass,
                     "job=i"    => \$nRastJobID,
                     "genome=s" => \$sGenomeName,
                     "igbdir=s" => \$sIGBdir);

# PRINT HELP
$sHelp = 1 unless ($sRastPass and $sRastUser and $nRastJobID and $sGenomeName and $sIGBdir);
if (!$res || $sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName -u <rastuser> -p <rastpass> -j <rastjob> -g <genomename> -i <igbdir>
   
   Arguments 
    -u --user <string>
      RAST server username
    -p --pass <string>
      RAST server password
    -j --job <integer>
      RAST job ID
    -g --genome <string>
      IGB genome name
    -i --igbdir <string>
      IGB root directory
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Check args
die "Error: IGB directory '$sIGBdir' does not exist\n" unless (-d $sIGBdir);
die "Error: incorrect job ID format\n" unless ($nRastJobID =~ /^\d+$/);

# Create new genome dir
my $sGenomeDir = "$sIGBdir/$sGenomeName";
mkdir($sGenomeDir) or die "Error: could not create genome dir: $!\n";

# Retrieve gff3 annotations and convert rRNA and tRNA annots
open GFFOUT, ">$sGenomeDir/$sGenomeName.gff3" or die "Error: can't open gff3 file '$sGenomeDir/$sGenomeName.gff3' for writing: $!\n";
open GFFIN, "$sSvrRetrieveJob $sRastUser $sRastPass $nRastJobID gff3 |" or die "Error: could not retrieve gff file for job '$nRastJobID'\n";
while (<GFFIN>){
   next if (/^\s*$/);
   if (/^ *#/){
      print GFFOUT $_;
      next;
   }
   my @asLine = split /\t/, $_, -1;
   $asLine[0] =~ s/\|/_/g;
   $asLine[2] = 'exon' if ($asLine[2] eq 'rRNA');
   $asLine[2] = 'exon' if ($asLine[2] eq 'tRNA');
   print GFFOUT join("\t", @asLine);
}
close GFFIN;
close GFFOUT;

# Retrieve fasta sequence
my ($flSeq, $sLocusID) = (0,"");
open FASTAOUT, ">$sGenomeDir/$sGenomeName.fasta" or die "Error: can't open fasta '$sGenomeDir/$sGenomeName.fasta' for writing: $!\n";
open FASTAIN, "$sSvrRetrieveJob $sRastUser $sRastPass $nRastJobID genbank |"
   or die "Error: could not retrieve sequence file for job '$nRastJobID'\n";
while (<FASTAIN>){
   next if (/^\s*$/);
   next if (/^ *#/);
   if (/^LOCUS(.*)\s+\d+ bp\s+DNA/){
      $sLocusID = $1;
      $sLocusID =~ s/^ +//;
      $sLocusID =~ s/\|/_/g;
      $sLocusID =~ s/\s+$//;
      print FASTAOUT ">$sLocusID\n";
   }
   if (/^ORIGIN/){
      $flSeq = 1;
      next;
   }
   if (/^\/\//){
      $flSeq = 0;
   }
   if ($flSeq){
      s/ //g;
      s/^\d+//;
      print FASTAOUT $_;
   }
}
close FASTAIN;
close FASTAOUT;

# Convert fasta sequence to twobit format
system("$sFaToTwoBit -noMask $sGenomeDir/$sGenomeName.fasta $sGenomeDir/$sGenomeName.2bit") == 0
   or die "Error: conversion of fasta to twobit format failed for job '$nRastJobID': $!\n";

# Convert the gff3 format to bed basic
open BEDOUT, ">$sGenomeDir/$sGenomeName.bed" or die "Error: can't open bed file '$sGenomeDir/$sGenomeName.bed' for writing: $!\n";
open GFF3, "$sGffToBed $sGenomeDir/$sGenomeName.gff3|" or die "Error: can't open gff3 file '$sGenomeDir/$sGenomeName.gff3':$!\n";
while (<GFF3>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\r\n]+$//;
   my @asLine = split /\t/, $_, -1;
   my $sAnnots = $asLine[3];
   
   # Parse ID string
   my %hAnnots;
   my @asAnnots = split /;/, $sAnnots;
   foreach my $sAnnot (@asAnnots){
      my ($key, $val) = split /=/, $sAnnot;
      $hAnnots{$key} = $val;
   }
   my $sID   = exists($hAnnots{ID}) ? $hAnnots{ID} : "ID not found";
   my $sName = exists($hAnnots{Name}) ? $hAnnots{Name} : "Name not found";
   
   # Append bed basic fields
   $asLine[3] = $sID;
   push @asLine, $sID;
   push @asLine, $sName;
   print BEDOUT join("\t", @asLine), "\r\n";
}
close GFF3;
close BEDOUT;

# Create the annots file
open ANNOTSOUT, ">$sGenomeDir/annots.xml" or die "Error: can't open annots.xml file '$sGenomeDir/annots.xml' for writing: $!\n";
print ANNOTSOUT "<files>\r\n";
print ANNOTSOUT "<file name=\"$sGenomeName.bed\" title=\"Annotation\" description=\"Gene annotations\" label_field=\"ID\" background=\"FFFFFF\" foreground=\"008000\" positive_strand_color=\"008000\" negative_strand_color=\"008000\" show2tracks=\"true\" direction_type=\"both\" max_depth=\"10\" name_size=\"12\" connected=\"true\" load_hint=\"Whole Sequence\"/>\r\n";
print ANNOTSOUT "</files>\r\n";
close ANNOTSOUT;

# Create the genome file
my @aaFastaLengths = get_fasta_lengths("$sGenomeDir/$sGenomeName.fasta");
open GENOMEOUT, ">$sGenomeDir/genome.txt" or die "Error: can't open genome.txt file '$sGenomeDir/genome.txt' for writing: $!\n";
foreach my $rContig (@aaFastaLengths){
   print GENOMEOUT join("\t", $rContig->[0], $rContig->[1]), "\r\n";
}
close GENOMEOUT;

# And finally, append the genome to the content.txt file
my %hContentIDs;
open CONTENTOUT, ">$sIGBdir/contents_new.txt" or die "Error: can't open '$sIGBdir/contents_new.txt' for writing: $!\n";
open CONTENT, "$sIGBdir/contents.txt" or die "Error: can't open '$sIGBdir/contents.txt': $!\n";
while (<CONTENT>){
   next if (/^\s*$/);
   next if (/^ *#/);
   print CONTENTOUT $_;
   my @asLine = split /\t/;
   $hContentIDs{$asLine[0]}++;
}
close CONTENT;
unless(exists $hContentIDs{$sGenomeName}){
   print CONTENTOUT "$sGenomeName\t$sGenomeName\r\n";
}
close CONTENTOUT;
system("mv -f $sIGBdir/contents.txt $sIGBdir/contents.bkp") == 0 or die "Error: can't replace contents.txt file for job '$nRastJobID'\n";
system("mv -f $sIGBdir/contents_new.txt $sIGBdir/contents.txt") == 0 or die "Error: can't replace contents.txt file for job '$nRastJobID'\n";


#################
## SUBROUTINES ##
#################

# get_fasta_lengths
#
# Return lengths of sequences in a multi-fasta file
sub get_fasta_lengths {
   my ($sInput) = @_;
   my @aaReturn;
   my $sFastaHeader = '';
   my $sFastaSeq    = '';
   open INPUT, "<$sInput" or die "Error: can't read the fasta file\n";
   while (<INPUT>){
      if (/^>/ or eof){
         if (eof){
            die "Error: file ends in fasta header without sequence\n" if (/^>/);
            $sFastaSeq .= $_;
         }
         if ($sFastaHeader){
            $sFastaHeader =~ s/^>//;
            $sFastaHeader =~ s/\s+$//;
            $sFastaHeader =~ s/[\n\r]+//g;
            $sFastaSeq    =~ s/\s//g;
            $sFastaSeq    =~ s/[\n\r]+//g;
            push @aaReturn, [($sFastaHeader, length($sFastaSeq))];
         }
         $sFastaHeader = $_;
         $sFastaSeq    = "";
      }
      else{
         next if (/^\s*$/);
         next if (/^ *#/);
         $sFastaSeq .= $_ if ($sFastaHeader);
      }
   }
   close INPUT;
   return(@aaReturn);
}
