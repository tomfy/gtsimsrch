#!/usr/bin/perl -w
use strict;

my %id1_id2d = ();
while(<>){
  next if(/^\s*#/);
  my @cols = split(" ", $_);
  my ($id1, $id2, $d) = @cols[0,1,5];
  if(!exists $id1_id2d{$id1}){
    $id1_id2d{$id1} = [$id2, $d];
  }
}

my @sorted_id1s = sort { $id1_id2d{$a}->[1] <=> $id1_id2d{$b}->[1] } keys %id1_id2d;

for my $id1 (@sorted_id1s){
  my $id2d = $id1_id2d{$id1};
  print "$id1  ", $id2d->[0], "  ", $id2d->[1], "\n";
}
