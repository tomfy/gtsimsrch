#!/usr/bin/perl -w
use strict;

my %id = ();
while(<>){
  my @cols = split(" ", $_);
  $id{$cols[0]}++;
  $id{$cols[2]}++;
  $id{$cols[3]}++;
}

print join("\n", keys %id), "\n";
