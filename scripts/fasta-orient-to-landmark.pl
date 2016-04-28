#!/usr/bin/perl

# 16.10.2014 08:37:06 EDT
# Harm van Bakel <hvbakel@gmail.com>
use lib("/sc/orga/projects/InfectiousDisease/tools/File-Which-1.21/lib/", "/sc/orga/projects/InfectiousDisease/tools/File-Temp-1.21/lib");
# MODULES
use strict;
use warnings;
use Getopt::Long;
use File::Which;
use File::Temp qw/ tempfile tempdir /;

# GLOBALS
$ENV{TMPDIR}  ||= "/tmp";

# GET PARAMETERS
my $sHelp            = 0;
my $sGenomeFile      = "";
my $sLandmarkFile    = "";
my $sLandmarkSeqType = "dna";
my $sHeaderKey       = "";
my $sSuffix          = "";
my $nFlankSize       = 0;
my $nMatchThresh     = 0.9;
my $flSplit          = 0;
GetOptions("help!"          => \$sHelp,
           "genome:s"       => \$sGenomeFile,
           "landmark:s"     => \$sLandmarkFile,
           "type:s"         => \$sLandmarkSeqType,
           "key:s"          => \$sHeaderKey,
           "orientsuffix:s" => \$sSuffix, 
           "flank:n"        => \$nFlankSize,
           "matchlength:n"  => \$nMatchThresh,
           "split!"         => \$flSplit);

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
      Fasta file containing the landmark sequence(s)
    -t --type <string>
      The landmark sequence type Type is one of:
                 dna - DNA sequence
                 rna - RNA sequence
                 prot - protein sequence
      default: $sLandmarkSeqType
    -k --key <string>
      Optional string that should be present in the fasta header
      of the genome sequence, e.g. "circ". Can be used to ensure that
      only circularized genomes are re-oriented.
    -o --orientsuffix <string>
      Optional suffix to add to the fasta header of a re-oriented sequence
    -m --matchlength <float>
      Match length threshold as a fraction of the landmark sequence length
      Default: $nMatchThresh;
    -f --flank <integer>
      Length of the flanking sequence beyond the landmark at which the
      breakpoint will be positioned. Default: $nFlankSize
    -s --split
      Split rather than reorient contig at landmark
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Check arguments
$sLandmarkSeqType = lc($sLandmarkSeqType);
die "Error: flank size must be an integer\n"        unless ($nFlankSize =~ /^\d+$/);
die "Error: file '$sGenomeFile' does not exist\n"   unless (-e $sGenomeFile);
die "Error: file '$sLandmarkFile' does not exist\n" unless (-e $sLandmarkFile);
die "Error: unknown landmark sequence type '$sLandmarkSeqType'. Use either 'dna', 'rna', or 'prot'.\n" unless ($sLandmarkSeqType =~ /^(dna|rna|prot)$/);

# Check prerequisites
my $sBlatBin       = which('blat');
my $sFaToTwoBitBin = which('faToTwoBit');
die "Error: can't find blat binary in path\n"     unless ($sBlatBin);
die "Error: can't find faToTwoBit binary in path" unless ($sFaToTwoBitBin);

# Get landmark best hit coordinates and orientation
my $sGenome2bitFile = fasta_to_twobit($sFaToTwoBitBin, $sGenomeFile);
my ($sLandmarkSeqID, $nLandmarkPos, $flRevComp) = get_landmark_position($sBlatBin, $sGenome2bitFile, $sLandmarkFile, $nFlankSize, $nMatchThresh);

# Make sure the suffix is prefixed with a pipe character
if ($sSuffix){
   $sSuffix =~ s/^\s*_*\|*//;
   $sSuffix = "|${sSuffix}";
}

# Process genomic fasta file
if ($sLandmarkSeqID){
   # We got a hit to the landmark sequence, proceed with reorientation
   my ($raGenomeSequences, $nFaLineSize) = read_fasta($sGenomeFile);
   for (my $i=0 ; $i<@$raGenomeSequences ; $i++){
      my ($sID, $sSeq) = @{$raGenomeSequences->[$i]};
      
      # Reorient sequence if we found the hit ID
      my $flProceed = 0;
      if ($sID eq $sLandmarkSeqID){
         # Check the header key if requested
         if ($sHeaderKey){
            if ($sID =~ /$sHeaderKey/){
               $flProceed = 1;
            }
         }
         else{
            $flProceed = 1;
         }
      }
            
      # Reorient contig if conditions are met
      if ($flProceed){
         $sSeq = reverse_complement($sSeq) if ($flRevComp);
         if ($flSplit){
            my $sSeqA = substr($sSeq, $nLandmarkPos);
            $sSeqA =~ s/.{$nFaLineSize}/$&\n/sg;
            $sSeqA =~ s/\n+$//;
            print ">${sID}${sSuffix}_splitA\n${sSeqA}\n";
            
            my $sSeqB = substr($sSeq, 0, $nLandmarkPos);
            $sSeqB =~ s/.{$nFaLineSize}/$&\n/sg;
            $sSeqB =~ s/\n+$//;
            print ">${sID}${sSuffix}_splitB\n${sSeqB}\n";
         }
         else{
            $sSeq = substr($sSeq, $nLandmarkPos) . substr($sSeq, 0, $nLandmarkPos);
            $sSeq =~ s/.{$nFaLineSize}/$&\n/sg;
            $sSeq =~ s/\n+$//;
            print ">${sID}${sSuffix}\n${sSeq}\n";
         }
         warn "Reoriented sequence '$sID' to landmark position $nLandmarkPos\n\n";
      }
      else{
         $sSeq =~ s/.{$nFaLineSize}/$&\n/sg;
         $sSeq =~ s/\n+$//;
         print ">$sID\n$sSeq\n";
      }
   }
}
else{
   # We didn't get a hit so we just output the original sequence
   warn("Warning: skipped reorientation of '$sGenomeFile' since it has no hit to landmark sequence.\n\n");
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
   
   if ($sLandmarkSeqType eq 'prot'){
      system("$sBlatBin -t=dnax -q=prot -minIdentity=90 -noHead $sGenome2bitFile $sLandmarkFile $sBlatOutput > /dev/null") == 0 or die "Error: blat search failed: $!\n";
   }
   else{
      system("$sBlatBin -t=dna -q=dna -minIdentity=90 -noHead $sGenome2bitFile $sLandmarkFile $sBlatOutput > /dev/null") == 0 or die "Error: blat search failed: $!\n";
   }
   
   # Read blat output
   my %hHitPositions;
   my ($nTopHitScore, $sTopHitSeqID, $nTopHitPos, $flTopRevComp) = (0, "", 0, 0);
   warn("BLAT HITS FOR LANDMARK(S)\n-------------------------\n");
   open BLAT, "$sBlatOutput" or die "Error: can't open blat output file: $!\n";
   while (<BLAT>){
      warn($_);
      next if (/^\s*$/);
      next if (/^ *#/);
      s/[\n\r]+$//;
      my @asLine = split /\t/;
      my $nQsize    = $asLine[10];
      my $nTsize    = $asLine[14];
      my $nHitScore = $asLine[0] - $asLine[1] - $asLine[4] - $asLine[6];
      my $flRevComp = $asLine[8] =~ /\+$/ ? 0 : 1;
      my $sSeqID    = $asLine[13];
      
      # Get start coordinate of hit position, taking into account that we're going
      # to reverse-complement the genome in case of a negative strand match
      my $nHitPos   = $flRevComp ? $asLine[16] + $nFlankSize : $asLine[15] - $nFlankSize;
      $nHitPos = 0       if ($nHitPos < 0);
      $nHitPos = $nTsize if ($nHitPos > $nTsize);
      $nHitPos = $nTsize - $nHitPos if ($flRevComp);
      if ($nHitScore/$nQsize >= $nMatchThresh){
         $hHitPositions{$nHitPos}++;
         if ($nHitScore > $nTopHitScore){
            ($nTopHitScore, $sTopHitSeqID, $nTopHitPos, $flTopRevComp) = ($nHitScore, $sSeqID, $nHitPos, $flRevComp);
         }
      }
   }
   close BLAT;
   unlink($sBlatOutput);
   unlink($sGenome2bitFile);
   warn("\n");
   
   # Check how many distinct hit locations we found and report if those hits are more than 10 nucleotides apart
   if ($sTopHitSeqID){
      my $nHitPosOffset  = 0;
      my @anHitPositions = sort {$a <=> $b} (keys %hHitPositions);
      if (@anHitPositions){
         $nHitPosOffset = $anHitPositions[$#anHitPositions] - $anHitPositions[0];
      }
      warn("Warning: found multiple hit locations for landmark sequence - picking best hit\n") if ($nHitPosOffset>10);
   }
   
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
