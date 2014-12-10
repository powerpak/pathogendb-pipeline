#!/use/local/bin/perl -w
#######  #######  #######  #######  
## FastaReader
#######  #######  #######  #######  
use strict;
use Symbol;
use IO::File;
package common_util::generic::FastaReader;

my $RND=2.0;
# my ($P_cur_seq,$P_cur_def);
# my ($cur_seq,$cur_def);
# my $cur_def;

##### ##### ##### ##### ##### ##### ##### 
##### new constructor.
##### ##### ##### ##### ##### ##### #####
sub new{	
 my $name = shift;
 my $class = ref($name) || $name;
 my $this = {};
 bless $this,$class;
 my $div =  shift; ## the demarcation used to show new records.
 if(!defined($div)){ $div=">";}
 $this->{div} = $div;
 $this->{cur_def}=undef;
 $this->{next_def} = undef;
 return $this;
}

sub uniquify_Fa{
    my $this=shift;
    my $P_files=shift;
    my $ofile=shift;
    if(!defined($ofile)){ 
	die "Usage:$0 file_list_reference ofile"; 
    } 
    print STDERR "combining @{$P_files}\ninto\n$ofile\n";
    open OUT,">$ofile" or die "ERROR: uniquify_file failed to open $ofile";
    my %acc_done;
    foreach my $file(@{$P_files}){ 
	$this->init_file($file);
	while(my ($P_defn,$P_body)=$this->next()){ 
	  my $seqname = $$P_defn;
	  $seqname =~ s/\s.*//;
	    if(defined($acc_done{$seqname})){next;}
	    $acc_done{$seqname}=1;
	    print OUT "$$P_defn\n";
	    print OUT "$$P_body\n";
	} 
	$this->close_file();
    } 
    close OUT;
    return;
} 

#### ################ ############# ###############
#### mergeFiles
### assumes that all the files are in same order etc.
### inefficient use of memory, gets everything in and 
### then does the merge.
#### can be changed, later,
#### ################ ############# ###############
sub mergeFiles{
    my $this=shift;
    my $P_files=shift;
    my $outfile=shift;
    my $fnum=@{$P_files};
    my @defns;my @bodys;
    for (my $i=0;$i<$fnum;$i++) {
	$this->init_file($P_files->[$i]);
	my $cnt=0;
	while (my ($P_defn,$P_body)=$this->raw_next()) {
	    if (!defined($defns[$cnt])){
		$$P_defn=~s/len=\d+\s*$//;
		$defns[$cnt]=$$P_defn;
	    }
	    elsif($defns[$cnt] !~ /$$P_defn/)
	      {die "wrong in ".$P_files->[$i]."<$$P_defn>ne<$defns[$cnt]>";}
	    if (!defined($bodys[$cnt])){$bodys[$cnt]=$$P_body;}
	    else {$bodys[$cnt].=$$P_body;}
	    $cnt++;
	}
	$this->close_file();
    }
    my $dnum=@defns;
    open OUT,">$outfile" or die "$outfile";
    for (my $i=0;$i<$dnum;$i++) {
	print OUT "$defns[$i]\n";
	print OUT "$bodys[$i]";
    }
    close OUT;
    return 1;
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### split_file_rand
#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
sub split_file_rand{
  my $this = shift;
  my $file=shift;
  my $P_ofiles=shift; ## [[file1,0.1],[file2,0.3] etc.] ##total of prob 
  ## should add upto 1, or else the last guy will get the diff

  if(!defined($file) || !defined($P_ofiles)){ 
    die "USage fr->split_file_rand(filename,[[ofile1,0.1],[ofile2,0.9]])\n".
	"\t10 percetn will go to ofile1\n".
	"\t90 percent will go to ofile2\n".
	"def lines become cleaner etc.";
  }
  my $clean =shift;
  if(!defined($clean)){$clean=0;}

  # if($unit_size==1){$clean=1;};

  $this->init_file($file);


  my @files;
  my @rnd;
  my $rnd=0;
  for(my $i=0;$i<@{$P_ofiles};$i++){ 
      #XXXXXXXXXXXXX
      $files[$i]= main::gensym();
#      push(@files,$_->[0]);
      print STDERR "creating file ".$P_ofiles->[$i]->[0]."\n";
      open $files[$i],">".$P_ofiles->[$i]->[0] or die "cannot open".$P_ofiles->[$i]->[0]."\n";
      $rnd+=$P_ofiles->[$i]->[1];
      push(@rnd,$rnd);
  } 
  if($rnd[-1]!= 1){$rnd[-1]=1;}

  while(my ($P_defn,$P_body)=$this->raw_next()){ 
    if($clean==1){$this->fixDefBod($P_defn,$P_body);}

    my $rnd=rand();
    for(my $j=0;$j<@rnd;$j++){
	if($rnd<$rnd[$j]){ 
	    print {$files[$j]} "$$P_defn\n"; 
	    print {$files[$j]} "$$P_body"; 
	    last;
	} 
    } 
  }
  for(my $i=0;$i<@{$P_ofiles};$i++){ 
      close $files[$i] or die "cannot close:".$P_ofiles->[$i]."\n";
  } 

  close $this->{file};
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### split_file
#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
sub split_file{
  my $this = shift;
  my $file=shift;
  my $unit_size=shift;

  if(!defined($file) || !defined($unit_size)){ 
    return my_die("USage fr->split_file(filename,size,[clean =0/1],[outdir]]) where \n".
     "\tsize is number of elelments in the fragmented files created \n".
     "\tclean will result in ambiguous codes getting fixed,capitalization,\n".
     "def lines become cleaner etc.");
  }
  my $clean =shift;
  my $outdir=shift;
  if(!defined($clean)){$clean=0;}

  $this->init_file($file);

  my $tmpfile = $file;
  my ($dir,$ext,$fname)=$this->getDirExt($tmpfile);

  if(defined($outdir)){$tmpfile=~"$outdir/$fname";}
  else{$outdir=$dir;$tmpfile=~ s/\.(fa.*)$//;}

  my $filecnt=0;
  my $fcnt=0;
  $fname="${tmpfile}_$filecnt.$ext";

  ## you want to use the definition, if only one fasta per file
  if($unit_size==1){$fcnt=$unit_size+10;}
  else {
    if(-e $fname){ 
      return my_die("split_file:$fname already exists");
    } 
    open OUT,">$fname" or return my_die("split_file:FAIL:>$fname");
  } 

  while(my ($P_defn,$P_body)=$this->raw_next()){ 
    $fcnt++;
    if($clean==1){$this->fixDefBod($P_defn,$P_body);}
    if($fcnt>$unit_size){ 
      if($fcnt > ($unit_size+3)){ } 
      else {close OUT or return my_die("split_file:FAIL:$fname");}
      $fcnt=1;$filecnt++;
      if($unit_size==1){ 
	my ($tdefn)=$$P_defn=~/^>(\S+)/;
	$fname = "$outdir/$tdefn.$ext";
      } else{$fname="${tmpfile}_$filecnt.$ext";}
      open OUT,">$fname" or return my_die("split_file:FAIL:$fname");
    }
    print OUT "$$P_defn\n"; 
    print OUT "$$P_body"; 
  }
  close OUT or return(my_die("FAIL: split_file:$fname"));
  close $this->{file};
  return (1,"good");
}





#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### split_file2
### does the split into a number of files, without knowing how many files
###  exist to begin with
#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
sub split_file2{
  my $this = shift;
  my $file=shift;
  my $file_num=shift;
  if(!defined($file) || !defined($file_num)){ 
    return(my_die("USage fr->split_file2(filename,num_of_files,[clean =0/1]) where \n".
      "\tnum_of_files is number of fragmented files created \n".
	"\tclean will get ambiguous codes fixed, capitalization,\n".
	  "def lines become cleaner etc."));
  } 
  my $clean =shift;
  if(!defined($clean)){ $clean=0;}
  # if($unit_size==1){$clean=1;};
  $this->init_file($file);
  my $tmpfile = $file;
  #  my ($dir,$ext) = $tmpfile =~ /^(.+)\/[^\/]+\.(fa.*)$/;
  ## to allow for .qual
  #  my ($dir,$ext) = $tmpfile =~ /^(.+)\/[^\/]+\.([^\.]*)$/;

  my ($dir,$ext,$fname)=$this->getDirExt($tmpfile);

#  $tmpfile=~ s/\.(fa.*)$//;
  $tmpfile=~ s/\.$ext$//;

  my @fname;
  my @fh;
  for(my $i=0;$i<$file_num;$i++){
    $fname[$i]="${tmpfile}_$i.$ext";

    if(-e $fname[$i]){ 
      return my_die("split_file2:$fname[$i] already exists");
    } 
    $fh[$i] = IO::File->new(">$fname[$i]") or return my_die("FAIL:split_file2:$fname[$i]");
    # open $fh[$i],">$fname[$i]" or die "$fname[$i]";
  } 
  my $fcnt=0;
  while(my ($P_defn,$P_body)=$this->raw_next()){ 
    if($clean==1){$this->fixDefBod($P_defn,$P_body);}
    print {$fh[$fcnt]} "$$P_defn\n"; 
    print {$fh[$fcnt]} "$$P_body"; 
    $fcnt++;$fcnt %= $file_num;
  } 
  for(my $i=0;$i<$file_num;$i++){
    close $fh[$i] or return(my_die("split_file2:FAIL:$fname[$i]"));
  } 
  close $this->{file};
  return (1,"good");
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### cnt_n_file
#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
sub cnt_in_file{
  my $this = shift;
  my $file=shift;
  
  my $lc=`grep "^>" $file | wc `;
  chomp($lc);
  return $lc;
} 
## number of zeros to pad the name to allow alphabetical sorting
sub make_pad{
  
  my $cnt=shift;
  my $dig=shift;

  my $pad="";
  if($cnt<10){ 
    $dig--;
  } elsif($cnt<100){ 
    $dig-=2;
  } elsif($cnt<1000){ 
    $dig-=3;
  } 
  for(my $i=0;$i<$dig;$i++){ 
    $pad.="0";
  }
  return $pad;
} 

sub my_die{
  my $mess=shift;
  return(0,$mess);
} 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### split_file_n
## make n files, with the order preserved (so putting back the 
## files in numerical order leads to same order as original)
#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
sub split_file_n{
  my $this = shift;
  my $file=shift;
  my $n=shift;

  my $clean =shift;
  my $outdir=shift;

  if(!defined($clean)){$clean=0;}

  if(!defined($file) || !defined($n)){ 
    return my_die("USage:fr->split_file_n(filename,n,[clean =0/1],[outdir]]) where \n".
     "\tn is number of fragmented files created \n".
     "\tclean will result in ambiguous codes in headers getting fixed,capitalization,\n".
     "def lines become cleaner etc.");
  }

  my $dig=0;
  if($n<=10){ $dig=1;}  elsif($n<=100){ $dig=2;}  elsif($n<=1000){ $dig=3;}
  else { return my_die("split_file_n:do you really want to split file to $n pieces ? I refuse !!!!\n");}


  my $lc=$this->cnt_in_file($file);

  my $num_per_file=int(($lc/$n) + 0.5); 

  $this->init_file($file);

  my $tmpfile = $file;

  my ($dir,$ext,$fname)=$this->getDirExt($tmpfile);

  if(defined($outdir)){$tmpfile=~"$outdir/$fname";}
  else{$outdir=$dir;$tmpfile=~ s/\.(fa.*)$//;}

  my $filecnt=0;

  my $fcnt=0;

  my $pad=make_pad($filecnt,$dig);

  $fname="${tmpfile}_$pad$filecnt.$ext";

  ## you want to use the definition, if only one fasta per file
  #  if($unit_size==1){$fcnt=$unit_size+10;}
  #  else {open OUT,">$fname" or die ">$fname";}

  if(-e $fname){   return my_die("split_file_n:$fname already exists");    } 

  open OUT,">$fname" or return(my_die("FAIL:split_file_n:>$fname"));

  while(my ($P_defn,$P_body)=$this->raw_next()){ 
    $fcnt++;
    if($clean==1){$this->fixDefBod($P_defn,$P_body);}
    if($fcnt>$num_per_file && $filecnt < ($n-1)){ 
      close OUT or return(my_die(" split_file_n:$fname"));
      $fcnt=1;$filecnt++;
      $pad=make_pad($filecnt,$dig);
      $fname="${tmpfile}_$pad$filecnt.$ext";
      open OUT,">$fname" or return(my_die(" split_file_n:$fname"));
    }
    print OUT "$$P_defn\n"; 
    print OUT "$$P_body"; 
  }
  close OUT or return(my_die("split_file_n:close:$fname"));
  close $this->{file};
  return (1,"good");
}
#### #### #### #### #### #### #### #### #### #### #### #### #### #### 

##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### getDirExt
##### ##### ##### ##### ##### ##### ##### ##### ##### 
sub getDirExt{
    my $this=shift;
    my $tmpfile=shift;
    my ($ext)=$tmpfile =~ /\.([^\.]*)$/;
    my ($fname) = $tmpfile =~ /([^\/]+)\.?[^\.]*$/;
    my ($dir) = $tmpfile =~ /^(.+)\/[^\/]+\.?[^\.]*$/;
    print STDERR  "getDirExt dir=$dir, ext=$ext, fname=$fname\n";
    return ($dir,$ext,$fname);
}
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### fixDefBod
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
sub fixDefBod{
  my $this = shift;
  my $P_defn=shift;
  my $P_body=shift;
  if($$P_defn =~ /\|(\w+\.\d+)\|/){$$P_defn = ">$1";}
  elsif($$P_defn =~ /^>(\S+)/){$$P_defn = ">$1";}
  $$P_body =~ tr/a-z/A-Z/;
  $$P_body =~ tr/RYMKSWHBVD/GCCGGACGGG/; ## get rid of SNPS
}
#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### get_defns
## returns a string, 
##       if complicated defline, then >$simple_Def=>$complicatedDef in a line
##       else just >$defline in a line
#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
sub get_defns{
  my $this = shift;
  my $file=shift;
  if(!defined($file)){
    die "USage fr->get_defns(filename)--lists all defn lines returns string\n";
  } 
  $this->init_file($file);
  my $tmpfile = $file;
  my ($dir,$ext) = $tmpfile =~ /^(.+)\/[^\/]+\.(fa.*)$/;
  my $report="";
  while(my ($P_defn,$P_body)=$this->raw_next()){ 
    if($$P_defn =~ /\|(\w+\.\d+)\|/){$report .= ">$1=$$P_defn\n";}
    else{$report .="$$P_defn\n";} 
  } 
  close $this->{file};
  return $report;
} 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### close_file, initialise the file.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
sub close_file{my $this = shift; close $this->{file};} 

#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### init_file, initialise the file.
####  Give filename as stdin if you want to handle stdin.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
sub init_file{
    my $this=shift;
    my $name = shift;
    die "no input file to init_file in FastaReader" if(!defined($name));
    
    # my $file2 = main::gensym();
    my $file2;

## or else just keep the STDIN stuff, so that you can use legitimate filenames.
## useful in splitting very large files.
    if ($name=~/^stdin$/i) { ## to allow for stdin inputs, for large file handling
	$file2=\*STDIN;
    }else {
       $file2 = IO::File->new("<$name") or die "cannot open $name";
       # open $file2,"<$name" or die "cannot open $name";	
    }
    return $this->init_fileptr($file2);
}



#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### init_file_rand, initialize to pick random samples. give also rndome number
####  Give filename as stdin if you want to handle stdin.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
sub init_file_rand{
    my $this=shift;
    my $name = shift;
    my $rand=shift;
    die "no input file to init_file_rand in FastaReader" if(!defined($name));
    die "no probability input  to init_file_rand in FastaReader" if(!defined($rand));
    
    # my $file2 = main::gensym();
    my $file2;
    $RND=$rand;

## or else just keep the STDIN stuff, so that you can use legitimate filenames.
## useful in splitting very large files.
    if ($name=~/^stdin$/i) { ## to allow for stdin inputs, for large file handling
	$file2=\*STDIN;
    }else {
       $file2 = IO::File->new("<$name") or die "cannot open $name";
       # open $file2,"<$name" or die "cannot open $name";	
    }
    return $this->init_fileptr($file2);
}
#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### init_fileptr, initialise the file. using a filepointer instead of a name.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
sub init_fileptr{
  my $this=shift;
  my $file = shift;
  $this->{file} = $file;

  my $div = $this->{div};

  while($this->{next_def} = <$file>){

    ## trying to pick the next $div, so this should onlybe operational the first time around
#    if($div!~/\>/ && $this->{next_def} !~ /^$div/){next;}
    if($this->{next_def} !~ /^$div/){next;}

    if($this->{next_def}=~/^\s*#/){ next;}
    if($this->{next_def} !~ /^\s*$/){last;}
    else {$this->{next_def} = undef;}
  } 
  if(!defined($this->{next_def})){return 1;}
  chomp($this->{next_def});
  $this->{next_def}=~s/\s+$//; ## remove craze stuff at end, like \r, that chomp does not remove 

  if($this->{next_def} !~ /^$div/){
    print STDERR "FR says:error in FASTA file\n";
    close $this->{file};
    return 1;
  }

  return 0;
}
##################    ##################    ##################    
##### next -- the next definition and sequence  
##################    ##################    ##################    
sub next{
    my $this = shift;
    my ($P_cur_def,$P_cur_seq) = $this->raw_next();
    if(!defined($P_cur_seq)){return;}

    $$P_cur_seq =~ s/\s+//g; ## to take out all the spaces inside

    return($P_cur_def,$P_cur_seq);
}
##################    ##################    ##################    
##### raw_next -- the next definition and sequence(without removing any of
######             the formatting.  
##################    ##################    ##################    
sub raw_next{
  my $this = shift;
  if(!defined($this->{next_def})){return;} 
  $this->{cur_def} = $this->{next_def};
  $this->{next_def} = undef;
  my $cur_seq = "";
  my $file = $this->{file};
  my $div = $this->{div};
  while(defined(my $tmp = <$file>)){
    if($tmp=~/^$div/){
      chomp($tmp);
      ## removes crazy \r type characters I found in haifan's files from geo/sra 
      $tmp=~s/\s+$//;
      $this->{next_def} = $tmp;last;            
    }
    $cur_seq .= $tmp;
  }
  if(!defined($this->{next_def})){ close $this->{file};}
  my $cur_def = $this->{cur_def};

  if($RND<1.99){
      my $tmprnd=rand();
      # print "tmprnd=$tmprnd\n";
      if($tmprnd > $RND){ return $this->raw_next();} 
  } 
  return(\$cur_def,\$cur_seq);
}
1;

=head1 NAME

FastaReader - reads in FASTA files.

=head1 SYNOPSIS

     my $fasta_reader = Util::FastaReader->new();

     $fasta_reader->init_file("stats.fasta");

     while(my ($defn,$P_seq) = $fasta_reader->next){  }

     next returns pointer to the sequence and the definition line 
     and returns nothing at EOF


=head1 AUTHOR

Ravi Sachidanandam, CSHL. ravi@cshl.org


=cut


