#!/usr/bin/perl -w
use strict;



my %progid_count = ();
my %progid_lines = ();
my %progid_paridpairs = ();
while(<>){
  next if(/^\s*#/);
  my @cols = split(" ", $_);
  my $progid = $cols[0];
  $progid_count{$progid}++;
  $progid_lines{$progid} .= "  " . $_;
}

while(my($id, $count) = each %progid_count){
  my $the_line = $progid_lines{$id};
  my @cols = split(" ", $the_line);
  my $paridpair = $cols[2] . " " . $cols[3];
  print "$id $count \n", $progid_lines{$id} if($count > 1);
}
