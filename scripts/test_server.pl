#!/usr/bin/env perl

use RASTserver;
use strict;
use Getopt::Long;

my $usage="Usage: $0 [--url server-url] [--test] [--verbose] username passwd format jobid [ jobid jobid ... ]\n";

my $verbose;
my $is_test;
my $url;
my $n;
if(!GetOptions('verbose'=> \$verbose,
	       'test'=> \$is_test,
	       'url=s' => \$url))
{
    die $usage;
}

@ARGV > 3 or die $usage;

my $username=shift;
my $password=shift;
my $format=shift;
my @jobs=@ARGV;

my $opts={};
if($url)
{
    $opts->{-server}=$url;
}
if($is_test)
{
    $opts->{-test}=1;
}

my $rast=new RASTserver($username, $password, $opts);
$n=0;
while(1){
    sleep(20);
my $res=$rast->status_of_RAST_job({-job => \@jobs});
    $n++;
for my $job(@jobs)
{
    my $status_hash=$res->{$job};
    my $status=$status_hash->{status};
    my $err_msg="(error message: $status_hash->{error_msg})" if $status eq 'error';
    print STDERR "status for job $job: $status $err_msg\n";
    if($status=~/complete/){
#	my $res=$rast->retrieve_RAST_job({-job=>$job, -format=>$format, -filehandle=>\*STDOUT});
#	if($res->{status} eq 'error')
#	{
#	    die "Error retrieving job output: $res->{error_msg}\n";
#	}
#	print $res->{contents};
	exit;
    }
#    if($verbose)
#    {
#	for my $vs(@{$status_hash->{verbose_status}})
#	{
#	    my ($s,$v)=@$vs;
#	    print "$s: $v\n";

#	}
#    }
    
}
}
