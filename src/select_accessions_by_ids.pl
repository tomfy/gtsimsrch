#!/usr/bin/perl/ -w
use strict;

#

my $idfile = shift;
my $gtfile = shift;
my $idcol = shift // 1;
$idcol--;

my %id = ();
open my $fhin, "<", $idfile;
while(<$fhin>){
  next if(/^\s*#/);
  my @cols = split(" ", $_);
  if(scalar @cols >= $idcol){
    $id{$cols[$idcol]}++;
  }
}
close $fhin;

open $fhin, "<", $gtfile;
while(<$fhin>){ # skip initial comment lines
  next if(/^\s*#/);
  print if(/^\s*MARKER/); # print 1st non-comment line (with marker ids)
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
