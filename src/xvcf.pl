#!/usr/bin/perl -w
use strict;

while(<>){
  next if(/^\s*#/);
  my @cols = split(" ", $_);
  @cols = @cols[9..$#cols];
  for my $elem (@cols){
    my ($gt, $ds, $gp) = split(':', $elem);
    my ($p0, $p1, $p2) = split(',', $gp);
    print "$gt  $ds  $p0 $p1 $p2\n";
  }
}
