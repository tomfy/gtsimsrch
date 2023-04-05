#!/usr/bin/perl/ -w
use strict;

my $idfile = shift;
my $gtfile = shift;

my %id = ();
open my $fhin, "<", $idfile;
while(<$fhin>){
  if(/^\s*(\S+)/){
    $id{$1}++;
  }
}
close $fhin;

open $fhin, "<", $gtfile;
while(<$fhin>){
  next if(/^\s*#/);
  print;
  last;
}

while(<$fhin>){
  next if(/^\s*#/);
  if(/^\s*(\S+)/){
    if(exists $id{$1}){
      print;
    }
  }
}
