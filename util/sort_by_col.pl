#!/usr/bin/perl -w
use strict;
# use List::Util qw(min max sum);

my $col = shift // 1;
$col--;

my @values = ();
while(<>){
  next if(/^\s*(#|Chr)/);
  my @cols = split(" ", $_);
  my $v = $cols[$col];
  next if($v eq '-nan'  or  $v eq 'nan');
  push @values, $v;
}

@values = sort {$a <=> $b} @values;

print join("\n", @values);
