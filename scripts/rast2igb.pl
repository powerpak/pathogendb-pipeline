#!/usr/bin/perl

# 08.01.2015 06:07:03 PST
# Harm van Bakel <hvbakel@gmail.com>
#
# Downloads annotations and sequence from a submitted RAST job and converts the output
# files into an IGB-compatible Quickload directory.
# 
# Alternatively, you can start from a GenBank file already downloaded from RAST or elsewhere
# (see the --from option)
#
# For our purposes, we strip all "|" pipe characters from contig IDs (replaced with "_")
# and create a .2bit file and BED annotation track for each genome.
#
# Full specifications for an IGB Quickload directory's conventions can be found at:
# https://wiki.transvar.org/display/igbman/Creating+QuickLoad+Sites
#
# Required BioPerl >1.6. (on Minerva, you should `module load bioperl` before running this)

# MODULES
use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use List::MoreUtils qw( minmax );
use File::Copy;
use File::Basename;

# GLOBALS
my $sSvrRetrieveJob = 'svr_retrieve_RAST_job';
my $sFaToTwoBit     = 'faToTwoBit';
my %hTrackColors    = (read_spans_bin_forward => "369EAD",
                       read_spans_bin_reverse => "C24642",
                       read_starts_in_bin_f => "7F6084",
                       read_starts_in_bin_r => "86B402",
                       read_terminates_in_bin_f => "A2D1CF",
                       read_terminates_in_bin_r => "C8B631",
                       read_start_clipped_in_bin_f => "6DBCEB",
                       read_start_clipped_in_bin_r => "52514E",
                       read_end_clipped_in_bin_f => "4F81BC",
                       read_end_clipped_in_bin_r => "A064A1",
                       deletions_in_read => "369EAD",
                       insertions_in_read => "C24642",
                       total_reads => "000000");

# If SAS_DIR is set in the environment, use perl to interpret the plbins directly (preserving the environment)
if ($ENV{'SAS_DIR'}) { $sSvrRetrieveJob = "perl $ENV{'SAS_DIR'}/plbin/svr_retrieve_RAST_job.pl"; }

# GET PARAMETERS
my $sHelp            = 0;
my $sRastUser        = '';
my $sRastPass        = '';
my $nRastJobID       = '';
my $sGenomeName      = '';
my $sIGBdir          = '';
my $sFromGenbankFile = '';
my $sBigWigDir       = '';
my $sQcDir           = '';
my $sBamFile         = '';
my $res = GetOptions("help!"       => \$sHelp,
                     "user=s"      => \$sRastUser,
                     "pass=s"      => \$sRastPass,
                     "job=i"       => \$nRastJobID,
                     "from=s"      => \$sFromGenbankFile,
                     "wigdir=s"    => \$sBigWigDir,
                     "qcdir=s"     => \$sQcDir,
                     "bam=s"       => \$sBamFile,
                     "genome=s"    => \$sGenomeName,
                     "igbdir=s"    => \$sIGBdir);

# PRINT HELP
$sHelp = 1 unless (($sFromGenbankFile or ($sRastPass and $sRastUser and $nRastJobID)) and $sGenomeName and $sIGBdir);
if (!$res || $sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: 
       $sScriptName -f <genbankfile> -g <genomename> -i <igbdir> [ -b <bigwigdir> -q <qcdir> -b <bamfile> ]
       $sScriptName -u <rastuser> -p <rastpass> -j <rastjob> -g <genomename> -i <igbdir> [ -b <bigwigdir> -q <qcdir> -b <bamfile> ]
   
   Arguments
       -f --from <string>
         GenBank file that will be used as input (overrides -u, -p, and -j)
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
       -w --wigdir <string>
         Optional directory containing bigwig track files to add to annots file
       -q --qcdir <string>
         Optional directory with assembly QC data
       -b --bamfile <string>
         Optional bam file with aligned reads
       -help
         This help message
   
HELP
}


##########
## MAIN ##
##########

# Check args
die "FATAL: IGB directory '$sIGBdir' does not exist\n" unless (-d $sIGBdir);
if ($sBigWigDir){
   die "FATAL: BigWig track directory '$sBigWigDir' does not exist\n" unless (-d $sBigWigDir);
}
if ($sQcDir){
   die "FATAL: Assembly QC directory '$sQcDir' does not exist\n" unless (-d $sQcDir);
}
if ($sBamFile){
   die "FATAL: Bam file '$sBamFile' does not exist\n" unless (-e $sBamFile);
   die "FATAL: Bam file '$sBamFile' does not have a matching index file\n" unless (-e "$sBamFile.bai");
}
die "FATAL: incorrect job ID format\n" unless ($nRastJobID =~ /^\d+$/) || $sFromGenbankFile;
$sGenomeName =~ s/ /_/g;

# Create a new Quickload genome dir
my $sGenomeDir = "$sIGBdir/$sGenomeName";
if (-d $sGenomeDir) {
   die "FATAL: $sGenomeDir already exists, possibly from a failed run.\n       Please move or delete it before re-running this script.\n";
}
mkdir($sGenomeDir) or die "FATAL: could not create $sGenomeDir - $!\n";


#--------------------------------------------#
# Create annotation file in beddetail format #
#--------------------------------------------#

# Create a beddetail file for the annotations
if ($sFromGenbankFile) {
   genbank_to_beddetail($sFromGenbankFile, $sGenomeDir, $sGenomeName);
} 
else {
   rast_to_beddetail($sSvrRetrieveJob, $sRastUser, $sRastPass, $nRastJobID, $sGenomeDir, $sGenomeName);
}


#--------------------------------------------------------------#
# Create the annots.xml file with annotation and bigwig tracks #
#--------------------------------------------------------------#

# Write annotations track to annots.xml file
open ANNOTSOUT, ">$sGenomeDir/annots.xml" or die "Error: can't open annots.xml file '$sGenomeDir/annots.xml' for writing: $!\n";
print ANNOTSOUT "<files>\r\n";
print ANNOTSOUT "   <file name=\"$sGenomeName.bed\" title=\"Annotation\" description=\"Gene annotations\" label_field=\"ID\" background=\"FFFFFF\" foreground=\"008000\" positive_strand_color=\"008000\" negative_strand_color=\"008000\" show2tracks=\"true\" direction_type=\"both\" max_depth=\"10\" name_size=\"12\" connected=\"true\" load_hint=\"Whole Sequence\"/>\r\n";

# Gather list of bigwig track files and create a bigwig track dir if any are found
my @asBigWigFiles = ();
if ($sBigWigDir){
   opendir my($dir), $sBigWigDir or die "Can't open $sBigWigDir : $!\n";
   @asBigWigFiles = grep { /^.*bw$/ } readdir $dir;
   mkdir("$sGenomeDir/bigwig") or die "FATAL: could not create $sGenomeDir/bigwig - $!\n";
}

# Copy the bigwig track files and add entries to annots file
foreach my $sBigWigFile (@asBigWigFiles){
   $sBigWigFile =~ s/.bw$//;
   copy("$sBigWigDir/$sBigWigFile.bw", "$sGenomeDir/bigwig/$sBigWigFile.bw") or die "FATAL: could not copy '$sBigWigDir/$sBigWigFile.bw' to '$sGenomeDir/bigwig/$sBigWigFile.bw' - $!\n";
   copy("$sBigWigDir/$sBigWigFile.bwt", "$sGenomeDir/bigwig/$sBigWigFile.bwt") or die "FATAL: could not copy '$sBigWigDir/$sBigWigFile.bwt' to '$sGenomeDir/bigwig/$sBigWigFile.bwt' - $!\n";
   my $sForeground = exists($hTrackColors{$sBigWigFile}) ? $hTrackColors{$sBigWigFile} : "008000";
   print ANNOTSOUT "   <file name=\"bigwig/$sBigWigFile.bw\" title=\"$sBigWigFile\" description=\"$sBigWigFile\" background=\"FFFFFF\" foreground=\"$sForeground\"/>\r\n";
}

# Add a bam file
if ($sBamFile){
   my $sBamBasename = basename($sBamFile);
   copy($sBamFile, "$sGenomeDir/$sBamBasename") or die "FATAL: could not copy '$sBigWigDir/$sBamFile' to '$sGenomeDir/$sBamBasename' - $!\n";
   copy("$sBamFile.bai", "$sGenomeDir/$sBamBasename.bai") or die "FATAL: could not copy '$sBigWigDir/$sBamFile.bai' to '$sGenomeDir/$sBamBasename.bai' - $!\n";  
   print ANNOTSOUT "   <file name=\"$sBamBasename\" title=\"$sBamBasename\" description=\"$sBamBasename\"/>\r\n";
}

# Close annots file
print ANNOTSOUT "</files>\r\n";
close ANNOTSOUT;


#----------------------------------------------------------------------------------------#
# Get genome fasta sequence, convert to 2bit format and create the IGB 'genome.txt' file #
#----------------------------------------------------------------------------------------#

# Retrieve the RAST genbank file and extract fasta sequences
my ($flSeq, $sLocusID) = (0,"");
open FASTAOUT, ">$sGenomeDir/$sGenomeName.fasta" or die "Error: can't open fasta '$sGenomeDir/$sGenomeName.fasta' for writing: $!\n";
open FASTAIN, ($sFromGenbankFile || "$sSvrRetrieveJob $sRastUser $sRastPass $nRastJobID genbank |")
   or die "Error: could not retrieve sequence file for job '$nRastJobID'\n";
while (<FASTAIN>){
   next if (/^\s*$/);
   next if (/^ *#/);
   if (/^LOCUS(.*)\s+\d+ bp\s+D?NA/){
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

# Convert fasta sequences to twobit format
system("$sFaToTwoBit -noMask $sGenomeDir/$sGenomeName.fasta $sGenomeDir/$sGenomeName.2bit") == 0
   or die "Error: conversion of fasta to twobit format failed for job '$nRastJobID': $!\n";

# Create the genome file
my @aaFastaLengths = get_fasta_lengths("$sGenomeDir/$sGenomeName.fasta");
open GENOMEOUT, ">$sGenomeDir/genome.txt" or die "Error: can't open genome.txt file '$sGenomeDir/genome.txt' for writing: $!\n";
foreach my $rContig (@aaFastaLengths){
   print GENOMEOUT join("\t", $rContig->[0], $rContig->[1]), "\r\n";
}
close GENOMEOUT;


#-------------------------------------------------------------#
# Add the IGB Quickload dir to the top-level content.txt file #
#-------------------------------------------------------------#

# Append the new IGB Quickload dir to the content.txt file
my %hContentIDs;
open CONTENTOUT, ">$sIGBdir/contents_new.txt" or die "Error: can't open '$sIGBdir/contents_new.txt' for writing: $!\n";
`touch $sIGBdir/contents.txt`;
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

#----------------------------------------------#
# Stash QC data in the IGB folder if available #
#----------------------------------------------#

if ($sQcDir){
   system("cp -R $sQcDir/* $sGenomeDir") == 0 or die "FATAL: Could not copy assembly QC data to '$sGenomeDir' - $!\n";
}

#---------------------------------------#
# Try to get MLST data for this genome  #
#---------------------------------------#

# Fetch MLST info
if ($sGenomeDir =~ /difficile/i){
   system("fetch_mslt.py --fasta $sGenomeDir/$sGenomeName.fasta --mlst cdifficile --output $sGenomeDir/mlst.txt") == 0 or die "Error: fetch MLST info for job '$nRastJobID'\n";
}
else{
   system("fetch_mslt.py --fasta $sGenomeDir/$sGenomeName.fasta --mlst mlst --output $sGenomeDir/mlst.txt")  == 0 or die "Error: fetch MLST info for job '$nRastJobID'\n";
}


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

# genbank_to_beddetail
#
# Convert genbank file to beddetail format
sub genbank_to_beddetail {
   my ($sFromGenbankFile, $sGenomeDir, $sGenomeName) = @_;

   my $inSeq = Bio::SeqIO->new(-file   => "<$sFromGenbankFile",
                               -format => 'Genbank');
   open BEDOUT, ">$sGenomeDir/$sGenomeName.bed" or die "Error: can't open bed file '$sGenomeDir/$sGenomeName.bed' for writing: $!\n";
   while (my $nextSeq = $inSeq->next_seq) {
      my $sLocusID = $nextSeq->display_id();
      $sLocusID =~ s/^ +//;
      $sLocusID =~ s/\|/_/g;
      $sLocusID =~ s/\s+$//;
      
      # Each CDS/RNA feature becomes a line in the BED file
      for my $feature ($nextSeq->get_SeqFeatures) {
         if ($feature->primary_tag =~ /^CDS|tRNA|rRNA$/i) {
            my $location = $feature->location;
            my $sStrand = $location->strand > 0 ? '+' : '-';
            my @anBlockStarts = ();
            my @anBlockSizes = ();
            my ($nStart, $nEnd, $sBlockStarts, $sBlockSizes, $sID, $sName, $sDescription);
            
            # Must convert GenBank locations into BED's chromStart, chromEnd, blockStarts and blockSizes.
            if ($location->isa("Bio::Location::SplitLocationI")) {
               ($nStart, $nEnd) = minmax(map { $_->start } $feature->location->sub_Location);
               $nStart = $nStart - 1;
               for my $subLocation ( sort {($a->start cmp $b->start) * $location->strand} $feature->location->sub_Location ) {
                  push @anBlockStarts, $location->strand > 0 ? $subLocation->start - 1 - $nStart : $nEnd - $subLocation->end;
                  push @anBlockSizes, $subLocation->end - $subLocation->start + 1;
               }
            } else {
               $nStart = $feature->location->start - 1;
               $nEnd = $feature->location->end;
               push @anBlockStarts, 0;
               push @anBlockSizes, $nEnd - $nStart;
            }
            
            # Figure out best names, IDs, and descriptions from feature tags.
            if ($feature->has_tag("db_xref")) {
               for my $sDbXref ($feature->get_tag_values('db_xref')) {
                  if ($sDbXref =~ /^SEED:(.+)$/) { $sID = $1; $sName = $1; }
               }
            }
            
            if ($feature->has_tag("gene")) {
               ($sID) = ($feature->get_tag_values('gene'));
            }
            if ($feature->has_tag("locus_tag")) {
               ($sName) = ($feature->get_tag_values('locus_tag'));
               ($sID)   = ($feature->get_tag_values('locus_tag')) unless ($sID);
            }
            
            die "Error: feature '$sID $sName' has no SEED ID\n" unless ($sID);
            if ($feature->has_tag("gene")) {
               $sName = ($feature->get_tag_values("gene"))[0];
            }
            if ($feature->has_tag("product")) {
               $sDescription = ($feature->get_tag_values("product"))[0];
            }
            if ($feature->primary_tag =~ /^tRNA/i) {
               $sName = $sDescription;
            } elsif ($feature->primary_tag =~ /^rRNA/i) {
               my @asDescParts = split /;\s+/, $sDescription;
               $sName = $asDescParts[-1];
            }
            
            # Print the line to the BED file.
            $sName =~ s/[^\w-]/-/g;
            print BEDOUT join("\t", ($sLocusID, $nStart, $nEnd, $sName, '0', $sStrand, $nStart, $nEnd, '0,0,0', 0+@anBlockSizes, 
                              join(',', @anBlockSizes), join(',', @anBlockStarts), $sID, $sDescription)) . "\n";
         }
      }
   }
   close BEDOUT;
}


# rast_to_beddetail
#
# Convert rast gff3 output to bedbasic format
sub rast_to_beddetail {
   my ($sSvrRetrieveJob, $sRastUser, $sRastPass, $nRastJobID, $sGenomeDir, $sGenomeName) = @_;

   # Get a gff3 formatted file from rast
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

   # Convert the GFF3 format to BED detail format
   my %hFeatures;
   open BEDOUT, ">$sGenomeDir/$sGenomeName.bed" or die "Error: can't open bed file '$sGenomeDir/$sGenomeName.bed' for writing: $!\n";
   open GFF, "sort -s -t '\t' -k9,9 $sGenomeDir/$sGenomeName.gff3 |" or die "Error sorting gff file: $!\n";
   while (<GFF>){
      next if /^\s*$/;
      next if /^\s*[Tt]rack/;
      next if /^\s*#/;
      s/[\n\r]$//g;
      my ($sChr, $sSource, $sType, $nStart, $nEnd, $nScore, $sStrand, $nFrame, $sID) = split /\t/;
      $sType = lc($sType);
      if (exists($hFeatures{$sID})){
         die "Error: feature '$sID' was found on multiple chromosomes\n" unless($hFeatures{$sID}{'chr'} eq $sChr);
         die "Error: feature '$sID' was found on multiple strands\n"     unless($hFeatures{$sID}{'strand'} eq $sStrand);
         $hFeatures{$sID}{$sType}{'start'} ||= $nStart;
         $hFeatures{$sID}{$sType}{'end'}   ||= $nEnd;
         $hFeatures{$sID}{$sType}{'start'} = $nStart if ($nStart < $hFeatures{$sID}{$sType}{'start'});
         $hFeatures{$sID}{$sType}{'end'}   = $nEnd   if ($nEnd   > $hFeatures{$sID}{$sType}{'end'});
         push @{$hFeatures{$sID}{$sType}{'features'}}, [$nStart, $nEnd];
      }
      else{
         print BEDOUT get_bedline(\%hFeatures) if (keys(%hFeatures));
         %hFeatures = ();
         $hFeatures{$sID}{'strand'}        = $sStrand;
         $hFeatures{$sID}{'chr'}           = $sChr;
         $hFeatures{$sID}{$sType}{'start'} = $nStart;
         $hFeatures{$sID}{$sType}{'end'}   = $nEnd;
         push @{$hFeatures{$sID}{$sType}{'features'}}, [$nStart, $nEnd];
      }   
   }
   print BEDOUT get_bedline(\%hFeatures) if (keys(%hFeatures));
   close BEDOUT;
}


# get_bedline
#
# Convert gff features to bed detail lines
sub get_bedline{
   my $rhFeatures = shift @_;
   my $sReturn = "";
   
   foreach my $sID (keys(%$rhFeatures)){
      # Parse ID string
      my %hAnnots;
      my @asAnnots = split /;/, $sID;
      foreach my $sAnnot (@asAnnots){
         my ($key, $val) = split /=/, $sAnnot;
         $hAnnots{$key} = $val;
      }
      my $sOutID   = exists($hAnnots{ID}) ? $hAnnots{ID} : "ID not found";
      my $sOutDesc = exists($hAnnots{Name}) ? $hAnnots{Name} : "Name not found";
   
      if (exists($rhFeatures->{$sID}{'exon'}) and exists($rhFeatures->{$sID}{'cds'})){
         my $nStart      = $rhFeatures->{$sID}{'exon'}{'start'} - 1;
         my $nThickStart = $rhFeatures->{$sID}{'cds'}{'start'} - 1;
         $sReturn .= join ("\t", $rhFeatures->{$sID}{'chr'}, $nStart, $rhFeatures->{$sID}{'exon'}{'end'},
                     $sOutID, '0', $rhFeatures->{$sID}{'strand'}, $nThickStart, $rhFeatures->{$sID}{'cds'}{'end'},
                     '0', get_blocks($rhFeatures->{$sID}{'exon'}{'start'}, $rhFeatures->{$sID}{'exon'}{'features'}), $sOutID, $sOutDesc);
      }
      elsif (exists($rhFeatures->{$sID}{'exon'})){
         my $nStart      = $rhFeatures->{$sID}{'exon'}{'start'} - 1;
         my $nThickStart = $nStart;
         $sReturn .= join ("\t", $rhFeatures->{$sID}{'chr'}, $nStart, $rhFeatures->{$sID}{'exon'}{'end'},
                     $sOutID, '0', $rhFeatures->{$sID}{'strand'}, $nThickStart, $rhFeatures->{$sID}{'exon'}{'end'},
                     '0', get_blocks($rhFeatures->{$sID}{'exon'}{'start'}, $rhFeatures->{$sID}{'exon'}{'features'}), $sOutID, $sOutDesc);
      }
      elsif (exists($rhFeatures->{$sID}{'cds'})){
         my $nStart      = $rhFeatures->{$sID}{'cds'}{'start'} - 1;
         my $nThickStart = $nStart;
         $sReturn .= join ("\t", $rhFeatures->{$sID}{'chr'}, $nStart, $rhFeatures->{$sID}{'cds'}{'end'},
                     $sOutID, '0', $rhFeatures->{$sID}{'strand'}, $nThickStart, $rhFeatures->{$sID}{'cds'}{'end'},
                     '0', get_blocks($rhFeatures->{$sID}{'cds'}{'start'}, $rhFeatures->{$sID}{'cds'}{'features'}), $sOutID, $sOutDesc);
      }
      else{
         die "Error: no exon or CDS information was found for feature '$sID'\n";
      }
      $sReturn .= "\r\n";
   }
   return($sReturn);
}


# get_blocks
#
# Get block count, sizes and starts
sub get_blocks {
   my ($nGenomeStart, $raBlocks) = @_;
   my @aaBlocks = @$raBlocks;
   
   if (scalar(@aaBlocks)==1){
      my ($nStart,$nEnd) = @{$aaBlocks[0]};
      my $blocksize = $nEnd-$nStart+1;
      return join("\t", 1,"$blocksize,","0,");
   }
   else{
      
      # Now sort by start then end and process each block
      @aaBlocks = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @aaBlocks;
      my @anBlockStarts;
      my @anBlockSizes;
      foreach my $rBlock (@aaBlocks){
         my ($nStart, $nEnd) = @{$rBlock};
         push @anBlockStarts, $nStart-$nGenomeStart;
         push @anBlockSizes,  $nEnd-$nStart+1;
      }
      my $sBlockStarts = join(',', @anBlockStarts);
      my $sBlockSizes  = join(',', @anBlockSizes);
      return join("\t", scalar(@aaBlocks), "$sBlockSizes,","$sBlockStarts,");
   }
}


