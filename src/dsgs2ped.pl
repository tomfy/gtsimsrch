#!/usr/bin/perl -w
use strict;

# read in a dosage matrix file
# write the equivalent .ped file and .map file
# alleles are all taken to be L and H, so
# dosage 0 -> L L, Dosage 1 -> L H, 2 -> H H

my @marker_ids = ();
while(<>){
  next if(/^\s*#/);
  @marker_ids = split(" ", $_);
  my $x = shift @marker_ids;
  die if($x ne 'MARKER');
  last;
}

while(my($i, $mid) = each @marker_ids){
  my $chrom = 100;
  my $position = -100;
  if($mid =~ /^S([0-9]+)_([0-9]+)/){
    $chrom = $1;
    $position = $2;
  }
  print STDERR "$chrom  " , $mid , " 0 " , $position, "\n";
}



while(<>){
  my @cols = split(" ", $_);
  my $accid = shift @cols;
  print "$accid $accid 0 0 0 -9  ";
  for my $dsg (@cols){
    if($dsg eq '0'){
      print "L L ";
    }elsif($dsg eq '1'){
      print "L H ";
    }elsif($dsg eq 'X'){
      print "0 0 ";
    }elsif($dsg eq '2'){
      print "H H ";
    }else{
      die "unexpected dosage: $dsg \n";
    }
  }
  print "\n";
}
