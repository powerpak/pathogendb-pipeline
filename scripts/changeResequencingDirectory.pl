my $file=shift;
my $dir1=shift;
my $dir2=shift;

open(FH, $file);
while (<FH>){
    if ($_ =~ /\{\{REFERENCE_DIRECTORY\}\}/){
		$_ =~ s/\{\{REFERENCE_DIRECTORY\}\}/$dir1\/$dir2/g;
    }
    print $_;
}
