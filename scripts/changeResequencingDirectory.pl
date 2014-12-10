my $file=shift;
my $dir1=shift;
my $dir2=shift;
open(FH, $file);
while(<FH>){
    if($_=~/circularized_sequence/){
	$_=~s/\/sc\/orga\/scratch\/attieo02\/VRE\/VRE_1/$dir1/g;
	$_=~s/circularized_sequence\/VRE_1/$dir2/g;
    }
    print $_;
}
