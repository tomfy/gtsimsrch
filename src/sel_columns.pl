#!/usr/bin/perl -w
use strict;

my $colstr = shift // '1'; # e.g. '1,5,6'

my @selcols = split(",", $colstr);
@selcols = map($_ - 1, @selcols);

while(<>){
  next if(/^\s*#/);
  my @cols = split(" ", $_);
  for my $acol (@selcols){
    print $cols[$acol], " ";
  }
  print "\n";
}
