#!/usr/bin/perl -w
use strict;

my %id1_mind = ();

while(<>){
  next if(/^\s*#/);
  my @cols = split(" ", $_);
  my ($id1, $d) = @cols[0,5];
  if(!exists $id1_mind{$id1}){
    $id1_mind{$id1} = $d;
  }
}

while(my ($id1, $d) = each %id1_mind){
  print "$id1 $d\n";
}
