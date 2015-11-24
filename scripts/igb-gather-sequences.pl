#!/usr/bin/perl

# 06.11.2015 09:51:13 EST
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp        = 0;
my $sIgbBaseDir  = '/sc/orga/projects/InfectiousDisease/igb/';
my $sInputList   = '';
my $sSuffix      = '';
GetOptions("help!"    => \$sHelp,
           "igb:s"    => \$sIgbBaseDir,
           "list:s"   => \$sInputList,
           "suffix:s" => \$sSuffix);

# PRINT HELP
$sHelp = 1 unless($sIgbBaseDir and $sInputList);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName
   
   Gather specific fasta sequences from a selection of IGB folders.
   
   Arguments:
    -list <string>
      Tab-delimited list of IGB folders and fasta sequence IDs. Format:
      genome-dir-name <tab> <fasta-header-id>
    -igb <string>
      IGB base folder
      default: $sIgbBaseDir;
    -suffix <string>
      Suffix to add to fasta filename
      default: <none>
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Check if IGB dir exists
$sIgbBaseDir =~ s/\/$//;
die "Error: '$sIgbBaseDir' does not exist or is not a directory\n" unless (-d $sIgbBaseDir);

# Read the input file list
die "Error: '$sInputList' does not exist\n" unless (-e $sInputList);
my $rhInputList = read_input_list($sInputList);

# Collect fasta sequences
foreach my $sIgbGenomeDir (keys %$rhInputList){
   if (-e "$sIgbBaseDir/$sIgbGenomeDir"){
      # Get the run ID and valid pathogenID from the IGB folder
      my ($sFolderName, $sRunID, $sAddID, $sIsolateID) = parse_ids_from_foldername("$sIgbBaseDir/$sIgbGenomeDir");
      
      # Read genome fasta file
      if (-e "$sIgbBaseDir/$sIgbGenomeDir/$sIgbGenomeDir.fasta"){
         my ($rhFastaSeq, $nFaLineSize) = read_fasta("$sIgbBaseDir/$sIgbGenomeDir/$sIgbGenomeDir.fasta");
         
         open OUT, ">${sIgbGenomeDir}${sSuffix}.fasta" or die "Error: can't write to '${sIgbGenomeDir}${sSuffix}.fasta': $!\n";
         foreach my $sSeqID (keys %{$rhInputList->{$sIgbGenomeDir}}){
            if (exists $rhFastaSeq->{$sSeqID}){
               my @asSequences = @{$rhFastaSeq->{$sSeqID}};
               foreach my $sSeq (@asSequences){
                  $sSeq =~ s/.{$nFaLineSize}/$&\n/sg;
                  $sSeq =~ s/\n+$//;
                  print OUT ">$sIsolateID|$sSeqID\n$sSeq\n";
               }
            }
         }
         close OUT;
      }
      else{
         warn("Fasta file '$sIgbBaseDir/$sIgbGenomeDir/$sIgbGenomeDir.fasta' does not exist, skipping\n");
      }
   }
   else{
      warn("Genome dir '$sIgbBaseDir/$sIgbGenomeDir' does not exist, skipping\n");
   }
}


#################
## SUBROUTINES ##
#################

# read_input_list
#
# Read the input list and return hash with fasta file selection
sub read_input_list {
   my ($sInputList) = @_;

   # Read IDs from file
   my %hIDs;
   if ($sInputList){
      open LIST, "<$sInputList" or die "Error: can't read ID file\n";
      while (<LIST>){
         next if (/^\s*$/);
         next if (/^ *#/);
         s/[\n\r]+$//;
         my ($sGenomeDir, $sFastaID) = split /\t/;
         $sFastaID =~ s/^>+//;
         $hIDs{$sGenomeDir}{$sFastaID}++ if ($sFastaID);
      }
      close LIST;
   }
   return \%hIDs;
}


# parse_ids_from_foldername
#
# Get the isolate, extract and SMRT portal ID from the IGB folder name
sub parse_ids_from_foldername {
   my ($sIgbDir) = @_;
   my ($sRunID, $sAddID, $sIsolateID) = ('', '', '');
   
   $sIgbDir =~ s/\/$//;
   my (@asFolderPath) = split /\//, $sIgbDir;
   my $sFolderName = pop @asFolderPath;
   
   my (@asFolderName) = split /\_/, $sIgbDir;
   if (@asFolderName >= 3){
      $sRunID     = pop @asFolderName;
      $sAddID     = pop @asFolderName;
      $sIsolateID = pop @asFolderName;
   }
   return($sFolderName, $sRunID, $sAddID, $sIsolateID);
}

# read_fasta
#
# Reads content of a multifasta file into hash of arrays
# Returns the hash and fasta line length
sub read_fasta {
   my ($sFASTA) = @_;
   
   my %hFasta;
   my $sFastaHeader = '';
   my $sFastaSeq    = '';
   my $nFaLineSize  = 0;
   open FASTA, "<$sFASTA" or die "Error: can't read the fasta file\n";
   while (<FASTA>){
      s/[\n\r]+$//;
      if (/^>/){
         die "Error: file ends in fasta header without sequence\n" if (eof);
         $sFastaSeq  =~ s/\s//g;
         push @{$hFasta{$sFastaHeader}}, $sFastaSeq if ($sFastaHeader);
         
         # Reset for the next sequence
         $sFastaHeader = $_;
         $sFastaHeader =~ s/\s*$//;
         $sFastaHeader =~ s/^>\s*//;
         $sFastaSeq    = "";
      }
      elsif (eof){
         $sFastaSeq .= $_;
         $sFastaSeq  =~ s/\s//g;
         push @{$hFasta{$sFastaHeader}}, $sFastaSeq if ($sFastaHeader);
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
   return (\%hFasta, $nFaLineSize);
}

