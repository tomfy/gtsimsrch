#!/usr/bin/perl -w
use strict;

my $n_pedigrees = shift // 100;

my @accids = ();
while(<>){
  next if(/^\s*#/);
  next if(/^\s*MARKER/);
  if(/^\s*(\S+)/){
    push @accids, $1;
  }
}

my $n_so_far =0;
while(1){
  my $idx1 = int(rand($n_pedigrees));
  my $idx2 = int(rand($n_pedigrees));
  my $idx3 = int(rand($n_pedigrees));
  next if($idx2 == $idx1  or  $idx3 == $idx1);
  print  $accids[$idx1], "   ", $accids[$idx2], "  ", $accids[$idx3], "\n";
  $n_so_far++;
  last if($n_so_far >= $n_pedigrees);
}
