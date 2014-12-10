#!/usr/bin/perl5.10.1
use strict;
use File::Basename;
use lib(dirname(__FILE__));
use common_util::generic::FastaReader;

my $repo_dir = dirname(dirname(__FILE__));
my $fr=common_util::generic::FastaReader->new();
my $file=shift @ARGV;
my %Length;
my @data;
my %beginning;
my %end;
my %rotated_start;
my %flag;
my $file1=$file;

$file1=~s/.fasta//;
$file1=~s/.fa//;

system("$repo_dir/vendor/MUMmer3.23/nucmer --maxmatch --nosimplify $file $file -p ${file1}_itself; "
    . "$repo_dir/vendor/MUMmer3.23/show-coords ${file1}_itself.delta > ${file1}_itself.coords");

$fr->init_file($file);
my $i=0;

while(my ($P_defn, $P_body)=$fr->next()){
  $i++;
  $$P_defn=~s/>//;
  $Length{$$P_defn}=length($$P_body);
}

open(FH, "${file1}_itself.coords") or die "Can't open ${file}_itself.coords: $!\n";
my $i=0;
my $length=0;

while(<FH>){
  @data=split;
  if(/^\s+\d+/&& $data[0]==1 && $data[1]<$Length{$data[11]} && $data[4]==$Length{$data[11]}){
    $beginning{$data[11]}=$data[0]; #+int(($data[1]-$data[0])/2);
    $end{$data[11]}=$data[3]; #+int(($data[1]-$data[0])/2);
    $rotated_start{$data[11]}=int(0.5*($end{$data[11]}-$beginning{$data[11]}))+1;
    # $flag{$data[11]}=1;
    print $data[0]." ".$data[1]." ".$data[3]." ".$data[4]."\n";
  }
}

my $fr1=common_util::generic::FastaReader->new();

$fr1->init_file($file);
my @contigs = keys %beginning;
open(FH1, ">${file1}_circularized.fasta");

while(my ($P_defn, $P_body)=$fr1->next()){
  for(my $i=0; $i<@contigs; $i++){
    if($$P_defn eq ">".$contigs[$i]){
      print FH1 ">".$contigs[$i]."_circ\n";
      my $circularized_dna=substr($$P_body, $beginning{$contigs[$i]}, $end{$contigs[$i]}-$beginning{$contigs[$i]}+1);
      # print length($circularized_dna)."\n";
      # print length(substr($circularized_dna, $rotated_start{$contigs[$i]}))."\n";
      # print length(substr($circularized_dna, 0, $rotated_start{$contigs[$i]}))."\n";
      print FH1 substr($circularized_dna, $rotated_start{$contigs[$i]}).substr($circularized_dna, 0, 
          $rotated_start{$contigs[$i]});
      # print FH1 substr($$P_body, $beginning{$contigs[$i]}, $end{$contigs[$i]}-$beginning{$contigs[$i]}+1)."\n";
      $flag{$$P_defn} = 1;
    }
  }
  if(!$flag{$$P_defn}){
    print FH1 $$P_defn."\n";
    print FH1  $$P_body."\n";
  }
}

close(FH1);
