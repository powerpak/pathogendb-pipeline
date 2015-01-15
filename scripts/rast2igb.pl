#!/usr/bin/perl

# 08.01.2015 06:07:03 PST
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GLOBALS
my $sSvrRetrieveJob = 'svr_retrieve_RAST_job';
my $sFaToTwoBit     = 'faToTwoBit';

# If SAS_DIR is set in the environment, use perl to interpret the plbins directly (preserving the environment)
if ($ENV{'SAS_DIR'}) { $sSvrRetrieveJob = "perl $ENV{'SAS_DIR'}/plbin/svr_retrieve_RAST_job.pl"; }

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
$sGenomeName =~ s/ /_/g;

# Create a new Quickload genome dir
my $sGenomeDir = "$sIGBdir/$sGenomeName";
mkdir($sGenomeDir) or die "Error: could not create genome dir: $!\n";

# Retrieve gff3 annotations and convert rRNA and tRNA annots to a standard feature type
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

# Convert the gff3 format to bed basic
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

# Create the annots file
open ANNOTSOUT, ">$sGenomeDir/annots.xml" or die "Error: can't open annots.xml file '$sGenomeDir/annots.xml' for writing: $!\n";
print ANNOTSOUT "<files>\r\n";
print ANNOTSOUT "<file name=\"$sGenomeName.bed\" title=\"Annotation\" description=\"Gene annotations\" label_field=\"ID\" background=\"FFFFFF\" foreground=\"008000\" positive_strand_color=\"008000\" negative_strand_color=\"008000\" show2tracks=\"true\" direction_type=\"both\" max_depth=\"10\" name_size=\"12\" connected=\"true\" load_hint=\"Whole Sequence\"/>\r\n";
print ANNOTSOUT "</files>\r\n";
close ANNOTSOUT;

# Retrieve the RAST genbank file and extract fasta sequences
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

# And finally, append the new IGB Quickload dir to the content.txt file
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
