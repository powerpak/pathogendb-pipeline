#!/usr/bin/perl

# 29.01.2015 13:15:27 EST
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use DBI;
use Getopt::Long;

# GLOBALS
my $sDbConf = "$ENV{HOME}/.my.cnf";  # MySQL conf file with db password

# GET PARAMETERS
my $sHelp    = 0;
my $sIgbDir  = "";
GetOptions("help!"   => \$sHelp,
           "input:s" => \$sIgbDir);

# PRINT HELP
$sHelp = 1 unless($sIgbDir);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName
   
   Arguments:
    -i <string>
      IGB genome directory 
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Check if the IGB dir exists and has a genome file
$sIgbDir =~ s/\/$//;
die "Error: '$sIgbDir' does not exist or is not a directory\n" unless (-d $sIgbDir);
die "Error: '$sIgbDir' does not contain a 'genome.txt' file\n" unless (-e "$sIgbDir/genome.txt");

# Get the run ID and valid pathogenID from the IGB folder
my ($sFolderName, $sRunID, $sAddID, $sIsolateID) = parse_ids_from_foldername($sIgbDir);

# Check the identifiers
my $flIDerror = 0;
$flIDerror = 1 unless ($sRunID =~ /^\d+$/);
$flIDerror = 1 unless ($sAddID =~ /^\d+[A-Z]+$/);
$flIDerror = 1 unless ($sIsolateID =~ /^[A-Z]{2}\d{5}$/);
if ($flIDerror){
   die <<FOLDER
   Error: IGB folder does not conform to standard naming guidelines.
   
   Please format the folder name according to the following example:
      <species>_ER00023_1A_020225
   
   where:
    - ER00023 is the isolate ID (two capitals followed by 5 numbers)
    - 1A      is the stock/extract ID (one or more numbers, followed by one or capitals)
    - 020225  is the SMRT portal run ID for this particular assembly

FOLDER
}
my $sExtractID = join('.', $sIsolateID, $sAddID);

# Get the genome stats
my ($nTotalSize, $nMaxContigLength, $nN50length, $nContigCount) = get_stats_from_genomefile("$sIgbDir/genome.txt");

#--------------------#
# UPLOAD STARTS HERE #
#--------------------#

# Check if database connection file exists
die "Error: can't find database connection details file '$sDbConf'\n" unless (-e $sDbConf);

# Open database connection
my $dbh = DBI->connect("DBI:mysql:vanbah01_pathogens"
                     . ";mysql_read_default_file=$ENV{HOME}/.my.cnf"
                     . ';mysql_read_default_group=vanbah01_pathogens',
                       undef, undef) or die "Error: could not establish database connection ($DBI::errstr)";

# Check database to make sure the extract ID exists
my $sSQL   = "SELECT E.extract_ID FROM tExtracts E WHERE E.extract_ID=\"$sExtractID\"";
my $oQuery = $dbh->prepare($sSQL);
my $nCount = $oQuery->execute();
$oQuery->finish();
if ($nCount eq '0E0'){
   $dbh->disconnect();
   die "Warning: extract ID '$sExtractID' was not found in the database, no data was uploaded\n";
}

# Save data into tAssemblies
$sSQL   = "REPLACE INTO tAssemblies (extract_ID, assembly_ID, assembly_data_link, contig_count, contig_N50, contig_maxlength, contig_sumlength ) VALUES( ?, ?, ?, ?, ?, ? ,?)";
$oQuery = $dbh->prepare($sSQL);
$nCount = $oQuery->execute($sExtractID, $sRunID, $sFolderName, $nContigCount, $nN50length, $nMaxContigLength, $nTotalSize);
if ($nCount){
   print "Loaded assembly '$sRunID' for extract '$sExtractID' into pathogenDB\n";
}
else{
   warn "Error: Could not load assembly '$sRunID' for extract '$sExtractID' into pathogenDB\n";
}

# Save data into tSequencing_runs
$sSQL   = "REPLACE INTO tSequencing_runs (extract_ID, sequence_run_ID, sequencing_platform, paired_end, run_data_link) VALUES( ?, ?, ?, ?, ?)";
$oQuery = $dbh->prepare($sSQL);
$nCount = $oQuery->execute($sExtractID, $sRunID, 'Pacbio', 'No', "http://smrtportal.hpc.mssm.edu:8080/smrtportal/#/View-Data/Details-of-Job/$sRunID");
if ($nCount){
   print "Loaded sequencing run '$sRunID' for extract '$sExtractID' into pathogenDB\n";
}
else{
   warn "Error: Could not load sequencing run '$sRunID' for extract '$sExtractID' into pathogenDB\n";
}

# Disconnect from database
$dbh->disconnect();


#################
## SUBROUTINES ##
#################

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


# get_stats_from_genomefile
#
# Get total contig size, max contig size and assembly N50 from the genome file
sub get_stats_from_genomefile {
   my ($sGenome) = @_;

   # Get max length and sum
   my @anContigLengths;
   my ($nSumContigLength, $nMaxContigLength, $nContigCount) = (0, 0, 0);
   open GENOME, $sGenome or die "Error: can't open '$sGenome': $!\n";
   while (<GENOME>){
      next if (/^\s*$/);
      next if (/^ *#/);
      s/[\n\r]+$//;
      my ($sContigID, $nContigLength) = split /\t/;
      die "Error: '$sGenome' contains a non-numeric contig length value on line $.\n" unless ($nContigLength =~ /^\d+$/);
      $nMaxContigLength = $nContigLength if ($nContigLength > $nMaxContigLength);
      push @anContigLengths, $nContigLength;
      $nSumContigLength += $nContigLength;
      $nContigCount++;
   }
   close GENOME;
   
   # Get N50
   my ($nN50, $nCumSum) = (0, 0);
   foreach my $nLength (sort {$b <=> $a} @anContigLengths){
      $nCumSum += $nLength;
      if ( ($nCumSum/$nSumContigLength) > 0.5){
         $nN50 = $nLength;
         last;
      }
   }
   
   return($nSumContigLength, $nMaxContigLength, $nN50, $nContigCount);
}

