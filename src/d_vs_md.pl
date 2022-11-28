#!/usr/bin/perl -w
use strict;

my $ptout_file = shift;

my %id_md = ();

open my $fh, "<", "$ptout_file";

while(<$fh>){
  next if(/^\s*#/);
  my @cols = split(" ", $_);
  $id_md{$cols[0]} = $cols[1];
}
close $fh;

open $fh, "<", "$ptout_file";
while(<$fh>){
  next if(/^\s*#/);
  my @cols = split(" ", $_);
  my ($id1, $id2, $id3, $d) = @cols[0,2,3,17];
  my $md_total = ($id_md{$id1} // 0);
  $md_total += ($id_md{$id2} // 0);
  $md_total += ($id_md{$id3} // 0);
  print "$id1 $id2 $id3 $md_total $d\n";
}
