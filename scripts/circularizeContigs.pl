#!/usr/bin/perl

# GENERAL MODULES
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use lib(dirname(__FILE__));
use common_util::generic::FastaReader;

# SET REPO DIR
$ENV{REPO_DIR} ||=  dirname(dirname(__FILE__));

# GET ARGUMENTS
my $sHelp           = 0;
my $sInputFile      = "";
my $nMaxOverhangLen = 12000;
my $nMaxOverhangPct = 0.05;
GetOptions("help!"          => \$sHelp,
           "input:s"        => \$sInputFile,
           "lenoverhang:i"  => \$nMaxOverhangLen,
           "pctoverhang:s"  => \$nMaxOverhangPct);
           

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
    -lenoverhang <integer>
      Maximum length of overhang allowed at each contig end.
      default: $nMaxOverhangLen
    -pctoverhang <number>
      Maximum overlap length as percentage of contig.
      default: $nMaxOverhangPct
    -help
      This help message
   
HELP
}

##########
## MAIN ##
##########

# Check arguments
die "Error: lenoverhang must be an integer >= 0\n"          unless ($nMaxOverhangLen =~ /^\d+$/);
die "Error: pctoverhang must be a number between 0 and 1\n" unless ($nMaxOverhangPct =~ /^\d\.?\d+$/);

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
system("$ENV{REPO_DIR}/vendor/MUMmer3.23/nucmer --maxmatch --nosimplify $sInputFile $sInputFile -p ${sOutPrefix}_itself 2> ${sOutPrefix}_nucmer.log;" .
       "$ENV{REPO_DIR}/vendor/MUMmer3.23/show-coords -T ${sOutPrefix}_itself.delta > ${sOutPrefix}_itself.coords;");


# Process overlaps from nucmer data and define new beginning and end coordinates of circularizable contigs
my %hStart;
my %hEnd;
my %hRotatedStart;
my %hOverhangSum;
open(FH, "${sOutPrefix}_itself.coords") or die "Can't open ${sOutPrefix}_itself.coords: $!\n";
while(<FH>){
   next unless /^\d+/;
   s/[\n\r]+$//;
   my ($nQstart, $nQend, $nTstart, $nTend, $nQlen, $nTlen, $nPctIdent, $sQname, $sTname) = split /\t/;
   if ( ($sQname eq $sTname) && ($nQstart != $nTstart) && ($nQend != $nTend) ){
      if (exists $hContigLengths{$sQname}){
         my $nQoverhang = $nQstart - 1;
         my $nToverhang = $hContigLengths{$sQname} - $nTend;
         my $nMaxOverhang = ($nMaxOverhangPct * $hContigLengths{$sQname}) < $nMaxOverhangLen ? $nMaxOverhangPct * $hContigLengths{$sQname} : $nMaxOverhangLen;
         if ($nQoverhang <= $nMaxOverhang && $nToverhang <= $nMaxOverhang){
            if (exists $hOverhangSum{$sQname} ){
               if ($hOverhangSum{$sQname} < ($nQoverhang + $nToverhang)){
                  $hStart{$sQname}        = $nQstart;
                  $hEnd{$sQname}          = $nTstart;
                  $hRotatedStart{$sQname} = int(0.5*($hEnd{$sQname} - $hStart{$sQname})) + 1;
                  $hOverhangSum{$sQname}  = $nQoverhang + $nToverhang;
                  warn("Updated end overlap :  contig: ${sQname}   length:$hContigLengths{$sQname}   start:${nQstart}-${nQend}   end:${nTstart}-${nTend}   left-overhang:${nQoverhang}   right-overhang:${nToverhang}\n");
               }
            }
            else{
               $hStart{$sQname}        = $nQstart;
               $hEnd{$sQname}          = $nTstart;
               $hRotatedStart{$sQname} = int(0.5*($hEnd{$sQname} - $hStart{$sQname})) + 1;
               $hOverhangSum{$sQname}  = $nQoverhang + $nToverhang;
               warn("Found end overlap   :  contig:${sQname}   length:$hContigLengths{$sQname}   start:${nQstart}-${nQend}   end:${nTstart}-${nTend}   left-overhang:${nQoverhang}   right-overhang:${nToverhang}\n");
            }
         }
      }
      else{
         die "Error: could not determine length for contig '$sQname'\n";
      }
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

