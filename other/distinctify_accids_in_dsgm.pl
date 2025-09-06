#!/usr/bin/perl -w
use strict;

# make all accession ids distinct in a dsgm file.

# read a dsgm file from stdin
# if an accession id occurs on more than one line
# the 2nd and later occurrences will have __2, __3, etc. appended to the id.
# every line is output to stdout, with no change
# except for this modification of accession ids
# usage example:
# distinctify_accids_in_dsgm.pl < IITA.dsgm > IITA_distinctids.dsgm

my %id_count = ();

while(<>){
  if(/^\s*#/){
    print;
  }elsif(/^MARKER/){
    print;
  }elsif(/^CHROMOSOME/){
    print;
  }else{
    process_one_line($_);
    last;
  }
}
while(<>){
 process_one_line($_);
}

while(my($accid, $count) = each %id_count){
  if($count > 1){
    print STDERR "acc id: accid  count: $count \n";
  }
}
    
sub process_one_line{
  my $line = shift;
   if($line =~ /^\s*(\S+)/){
    my $accid = $1;
    my $count = $id_count{$accid} // 0; # number of previous occurrences of this id
    $count++;
    $id_count{$accid} = $count;
    if($count > 1){
      $accid .= "__" . $count;
      $line =~ s/^\s*\S+/$accid/;
    }
    print $line;
  }
}
