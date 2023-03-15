#!/usr/bin/perl -w
use strict;

while(<>){
  next if(/^\s*#/);
  my @cols = split(" ", $_);
  for my $i (9..$#cols){
    my @fields = split(':', $cols[$i]);
    for my $f (@fields){
      my @sfs = split(',', $f);
      for my $sfs (@sfs){
	print "$sfs  ";
      }
    }
    print "\n";
  }
}
