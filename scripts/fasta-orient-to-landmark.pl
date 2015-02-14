#!/usr/bin/perl

# 16.10.2014 08:37:06 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;
use File::Which;
use File::Temp qw/ tempfile tempdir /;

# GLOBALS
$ENV{TMPDIR}  ||= "/tmp";

# GET PARAMETERS
my $sHelp         = 0;
my $sGenomeFile   = "";
my $sLandmarkFile = "";
my $sHeaderKey    = "";
my $nFlankSize    = 0;
my $nMatchThresh  = 0.95;
GetOptions("help!"         => \$sHelp,
           "genome:s"      => \$sGenomeFile,
           "landmark:s"    => \$sLandmarkFile,
           "key:s"         => \$sHeaderKey,
           "flank:n"       => \$nFlankSize,
           "matchlength:n" => \$nMatchThresh);

# PRINT HELP
$sHelp = 1 unless($sGenomeFile and $sLandmarkFile);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName -g <genome> -l <landmark>
   
   Reorient and reverse-complement a circularized genome to
   a given landmark nucleotide sequence.
   
   Arguments:
    -g --genome <string>
      (Multi) fasta file containing the genome sequence to re-orient.
    -l --landmark <string>
      Fasta file containing the landmark sequence
    -k --key <string>
      Optional string that should be present in the fasta header
      of the genome sequence, e.g. "circ". Can be used to ensure that
      only circularized genomes are re-oriented.
    -m --matchlength <float>
      Match length threshold as a fraction of the landmark sequence length
      Default: $nMatchThresh;
    -f --flank <integer>
      Length of the flanking sequence beyond the landmark at which the
      breakpoint will be positioned. Default: $nFlankSize
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Check arguments
die "Error: flank size must be an integer\n"        unless ($nFlankSize =~ /^\d+$/);
die "Error: file '$sGenomeFile' does not exist\n"   unless (-e $sGenomeFile);
die "Error: file '$sLandmarkFile' does not exist\n" unless (-e $sLandmarkFile);

# Check prerequisites
my $sBlatBin       = which('blat');
my $sFaToTwoBitBin = which('faToTwoBit');
die "Error: can't find blat binary in path\n"     unless ($sBlatBin);
die "Error: can't find faToTwoBit binary in path" unless ($sFaToTwoBitBin);

# Get landmark best hit coordinates and orientation
my $sGenome2bitFile = fasta_to_twobit($sFaToTwoBitBin, $sGenomeFile);
my ($sLandmarkSeqID, $nLandmarkPos, $flRevComp) = get_landmark_position($sBlatBin, $sGenome2bitFile, $sLandmarkFile, $nFlankSize, $nMatchThresh);

# Process genomic fasta file
if ($sLandmarkSeqID){
   # We got a hit to the landmark sequence, proceed with reorientation
   my ($raGenomeSequences, $nFaLineSize) = read_fasta($sGenomeFile);
   for (my $i=0 ; $i<@$raGenomeSequences ; $i++){
      my ($sID, $sSeq) = @{$raGenomeSequences->[$i]};
      
      # Reorient sequence if we found the hit ID
      if ($sID eq $sLandmarkSeqID){
         if ($sHeaderKey){
            if ($sID =~ /$sHeaderKey/){
               $sSeq = reverse_complement($sSeq) if ($flRevComp);
               $sSeq = substr($sSeq, $nLandmarkPos) . substr($sSeq, 0, $nLandmarkPos);
               warn "Reoriented sequence '$sID' to landmark position $nLandmarkPos\n";
            }
            else{
               warn("Warning: skipped reorientation of sequence '$sID' since it did not contain the header key '$sHeaderKey'\n");
            }
         }
         else{
            $sSeq = reverse_complement($sSeq) if ($flRevComp);
            $sSeq = substr($sSeq, $nLandmarkPos) . substr($sSeq, 0, $nLandmarkPos);
            warn "Reoriented sequence '$sID' to landmark position $nLandmarkPos\n";
         }
      }
      $sSeq =~ s/.{$nFaLineSize}/$&\n/sg;
      $sSeq =~ s/\n+$//;
      print ">$sID\n$sSeq\n";
   }
}
else{
   # We didn't get a hit so we just output the original sequence
   warn("Warning: skipped reorientation of '$sGenomeFile' since it has no hit to landmark sequence.\n");
   open FILE, $sGenomeFile or die "Error: can't open '$sGenomeFile'\n";
   while (<FILE>){
      print;
   }
   close FILE;
}

#################
## SUBROUTINES ##
#################

# fasta_to_twobit
#
# Convert fasta file to twobit format
sub fasta_to_twobit {
   my ($sBin, $sFastaFile) = @_;
   my (undef, $sTwoBitFile) = tempfile("twobitXXXXXXXX", OPEN=>0, DIR=>$ENV{TMPDIR});
   system("$sBin $sFastaFile $sTwoBitFile.2bit") == 0 or die "Error: 2bit conversion failed: $!\n";
   return "$sTwoBitFile.2bit";
}

# get_landmark_position
#
# Returns the start location of the landmark in the genome
sub get_landmark_position {
   my ($sBlatBin, $sGenome2bitFile, $sLandmarkFile, $nFlankSize, $nMatchThresh) = @_;
   my (undef, $sBlatOutput) = tempfile("blatoutXXXXXXXX", OPEN=>0, DIR=>$ENV{TMPDIR});
   system("$sBlatBin -t=dna -q=dna -minIdentity=90 -noHead $sGenome2bitFile $sLandmarkFile $sBlatOutput > /dev/null") == 0 or die "Error: blat search failed: $!\n";
   
   # Read blat output
   my ($nTopHitScore, $sTopHitSeqID, $nTopHitPos, $flTopRevComp, $nHitCount) = (0, "", 0, 0, 0);
   open BLAT, "$sBlatOutput" or die "Error: can't open blat output file: $!\n";
   while (<BLAT>){
      next if (/^\s*$/);
      next if (/^ *#/);
      s/[\n\r]+$//;
      my @asLine = split /\t/;
      my $nQsize    = $asLine[10];
      my $nTsize    = $asLine[14];
      my $nHitScore = $asLine[0] - $asLine[1] - $asLine[4] - $asLine[6];
      my $flRevComp = $asLine[8] eq '+' ? 0 : 1;
      my $sSeqID    = $asLine[13];
      
      # Get start coordinate of hit position, taking into account that we're going
      # to reverse-complement the genome in case of a negative strand match
      my $nHitPos   = $flRevComp ? $asLine[16] + $nFlankSize : $asLine[15] - $nFlankSize;
      $nHitPos = 0       if ($nHitPos < 0);
      $nHitPos = $nTsize if ($nHitPos > $nTsize);
      $nHitPos = $nTsize - $nHitPos if ($flRevComp);
      if ($nHitScore/$nQsize >= $nMatchThresh){
         if ($nHitScore > $nTopHitScore){
            ($nTopHitScore, $sTopHitSeqID, $nTopHitPos, $flTopRevComp) = ($nHitScore, $sSeqID, $nHitPos, $flRevComp);
            $nHitCount++;
         }
      }
   }
   close BLAT;
   unlink($sBlatOutput);
   unlink($sGenome2bitFile);
   warn("Warning: found multiple hits for landmark sequence - picking best hit\n") if ($nHitCount>1);
   return($sTopHitSeqID, $nTopHitPos, $flTopRevComp);
}

# read_fasta
#
# Reads content of a multifasta file into an array
# Returns the array and fasta line length
sub read_fasta {
   my ($sFASTA) = @_;
   
   my @aaFasta;
   my $sFastaHeader = '';
   my $sFastaSeq    = '';
   my $nFaLineSize  = 0;
   open FASTA, "<$sFASTA" or die "Error: can't read the fasta file\n";
   while (<FASTA>){
      s/[\n\r]+$//;
      if (/^>/){
         die "Error: file ends in fasta header without sequence\n" if (eof);
         $sFastaSeq  =~ s/\s//g;
         push @aaFasta, ([$sFastaHeader,$sFastaSeq]) if ($sFastaHeader);
         
         # Reset for the next sequence
         $sFastaHeader = $_;
         $sFastaHeader =~ s/\s*$//;
         $sFastaHeader =~ s/^>\s*//;
         $sFastaSeq    = "";
      }
      elsif (eof){
         $sFastaSeq .= $_;
         $sFastaSeq  =~ s/\s//g;
         push @aaFasta, ([$sFastaHeader,$sFastaSeq]) if ($sFastaHeader);
         $nFaLineSize = length($_) if (length($_) > $nFaLineSize);
      }
      else{
         next if (/^\s*$/);
         next if (/^ *#/);
         $sFastaSeq .= $_ if ($sFastaHeader);
         $nFaLineSize = length($_) if (length($_) > $nFaLineSize);
      }
   }
   close FASTA;
   $nFaLineSize = 100 if ($nFaLineSize > 1000); # Make sure we return a sane linesize in case the original file wasn't wrapped
   return (\@aaFasta, $nFaLineSize);
}

# rev_comp
#
# Returns the reverse-complement of the supplied sequence
sub reverse_complement{
   my $seq = shift(@_);
   my $rev = reverse $seq;
   $rev =~ tr/ACGTacgt/TGCAtgca/;
   return $rev;
}
