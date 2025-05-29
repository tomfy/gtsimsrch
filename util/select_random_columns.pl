#!/usr/bin/perl -w
use strict;

# keep first col, and comment lines
# and a set of columns which
# are a fraction $keep_prob of all columns (approx.)

my $keep_prob = shift // 0.5;

my @keep_indices = ();
while (<>) {
  if (/^\s*#/) {
    print;
  } else {
    my @cols = split(" ", $_);
    while (my($i, $f) = each @cols) {
      if ($i == 0  or  rand() <= $keep_prob) {
	push @keep_indices, $i;
	print "$f ";
      }
    }
    print "\n";
    last;
  }
}

while (<>) {
  if (/^\s*#/) {
    print;
  } else {
    my @cols = split(" ", $_);
    for my $i (@keep_indices) {
      print $cols[$i] . " ";
    }
    print "\n";
  }
}
