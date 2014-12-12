#!/usr/bin/perl5.10.1

# GENERAL MODULES
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

# SET REPO AND LIB DIR
$ENV{REPO_DIR} ||= dirname(dirname(__FILE__));
use lib("$ENV{REPO_DIR}/scripts");
use common_util::generic::FastaReader;

# GET ARGUMENTS
my $sHelp        = 0;
my $sInputFile   = "";
GetOptions("help!"   => \$sHelp,
           "input:s" => \$sInputFile);

# PRINT HELP MESSAGE
$sHelp = 1 unless($sInputFile);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName -i <assembly.fasta>
   
   When given a fasta-formatted assembly file, this script will attempt
   to circularize each contig by examining contig ends for overlaps.
   Returns a multi-fasta file with circularized contigs. Contigs that
   cannot be circularized are left untouched.
   
   Arguments:
    -input <string>
      Fasta input file with assembled contigs. 
    -help
      This help message
   
HELP
}

##########
## MAIN ##
##########

# Derive prefix for output file
my $sOutPrefix =  $sInputFile;
$sOutPrefix    =~ s/.fasta$//;
$sOutPrefix    =~ s/.fa$//;

# Get contig lengths from input file
my %hContigLengths;
my $fr=common_util::generic::FastaReader->new();
$fr->init_file($sInputFile);
while(my ($P_defn, $P_body)=$fr->next()){
   $$P_defn =~ s/>//;
   $hContigLengths{$$P_defn}=length($$P_body);
}

# Run nucmer to identify overlap regions
die "Error: could not find nucmer binary in '$ENV{REPO_DIR}/vendor/MUMmer3.23/'. Check repository path\n" unless (-e "$ENV{REPO_DIR}/vendor/MUMmer3.23/nucmer");
die "Error: could not find show-coords binary in '$ENV{REPO_DIR}/vendor/MUMmer3.23/'. Check repository path\n" unless (-e "$ENV{REPO_DIR}/vendor/MUMmer3.23/show-coords");
system("$ENV{REPO_DIR}/vendor/MUMmer3.23/nucmer --maxmatch --nosimplify $sInputFile $sInputFile -p ${sOutPrefix}_itself;" .
       "$ENV{REPO_DIR}/vendor/MUMmer3.23/show-coords ${sOutPrefix}_itself.delta > ${sOutPrefix}_itself.coords;");

# Process overlaps from nucmer data and define new beginning and end coordinates of circularizable contigs
my %hStart;
my %hEnd;
my %hRotatedStart;
open(FH, "${sOutPrefix}_itself.coords") or die "Can't open ${sOutPrefix}_itself.coords: $!\n";
while(<FH>){
   my @aData = split;
   if( /^\s+\d+/ && ($aData[0]==1) && ($aData[1]<$hContigLengths{$aData[11]}) && ($aData[4]==$hContigLengths{$aData[11]}) ){
      $hStart{$aData[11]}        = $aData[0];
      $hEnd{$aData[11]}          = $aData[3];
      $hRotatedStart{$aData[11]} = int(0.5*($hEnd{$aData[11]} - $hStart{$aData[11]})) + 1;
      print join(" ", $aData[0], $aData[1], $aData[3], $aData[4], $aData[11], $hStart{$aData[11]}, $hEnd{$aData[11]}, $hRotatedStart{$aData[11]}), "\n";
   }
}
close(FH);

# Trim overlap at end of circularizable contigs, rotate genome by 180 degrees to set new start coordinate, and output results
open(FH1, ">${sOutPrefix}_circularized.fasta") or die "Can't open ${sOutPrefix}_circularized.fasta: $!\n";;
my $fr1 = common_util::generic::FastaReader->new();
$fr1->init_file($sInputFile);
while(my ($P_defn, $P_body) = $fr1->next() ){
   $$P_defn =~ s/>//;
   if (exists $hStart{$$P_defn}){
      my $sCircSeq = substr($$P_body, $hStart{$$P_defn}, $hEnd{$$P_defn}-$hStart{$$P_defn}+1);
      $$P_body     = substr($sCircSeq, $hRotatedStart{$$P_defn}).substr($sCircSeq, 0, $hRotatedStart{$$P_defn});
      $$P_defn    .= "_circ";
   }
   $$P_body =~ s/.{60}/$&\n/sg;  # This chunks the fasta output into 60 chars per line
   $$P_body =~ s/\n+$//;
   $$P_body .= "\n";
   print FH1 ">$$P_defn\n";
   print FH1 $$P_body;
}
close(FH1);
