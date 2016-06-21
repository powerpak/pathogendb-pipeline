use lib("/sc/orga/work/attieo02/Getopt-ArgParse-1.0.6/blib/lib/", "/sc/orga/work/attieo02/Moo-2.000002/blib/lib/", "/sc/orga/work/attieo02/Module-Runtime-0.014/blib/lib/", "/sc/orga/work/attieo02/Devel-GlobalDestruction-0.13/blib/lib/", "/sc/orga/work/attieo02/Sub-Exporter-Progressive-0.001011/blib/lib/");
use Getopt::ArgParse;
$ap=Getopt::ArgParse->new_parser(
    prog=>'extractCoordsFromMAF.pl',
    description=>'Extracts all start and stop coordinates from MAF file',
    epilog=>' ');
$ap->add_arg('--file', '-f', required=>1, help=>'Input file');
$ap->add_arg('--number', '-n', required=>1, help=>'Number of strains');
$args=$ap->parse_args();
$file=$args->file;
#print $file."\n";
$number=$args->number;
#print $number."\n";
#my $file=shift @ARGV;
#my $number=shift @ARGV;
open(FH, $file) or die "Can't open $file:$!\n";
my $i=0;
my $mult_flag;
while(<FH>){
    if(/^a/){
	my @data1=split;
	$data1[3]=~s/mult=//;
	if($data1[3] eq $number){
	    $i++;
	    $mult_flag=1;
	}else{
	    $mult_flag=0;
	}
    }
    if(/^s/&&$mult_flag){
	my @data=split;
	my @data1=split(/\./,$data[1]);
	$length{$data1[0]}+=$data[3];
#	$start{$data[1]}{$i}=$data[2];
#	$length{$data[1]}{$i}=$data[3];
    }
}
close(FH);
#print $i."\n";
foreach $data(keys %length){
    print $data." ".$length{$data}."\n";
}
#	for(my $j=0; $j<$i; $j++){
#	print $data." ".$j." ".$start{$data}{$j}." ".$length{$data}{$j}."\n";
#    }
#}
