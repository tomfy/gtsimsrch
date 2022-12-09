#!/usr/bin/perl -w
use strict;

my $d0 = shift // 0.03;
my $ratio = shift // 2.0;
my $x1 = shift // 0.018;
my $x2 = shift // 0.09;
my ($n1, $n2, $n3, $n4, $n5) = (0, 0, 0, 0, 0);

while (<>) {
  my @cols = split(" ", $_);
  if(scalar @cols > 35){
  my $best_d = $cols[19];
  my $nextbest_d = $cols[42];
  if ($best_d > $d0) {
    $n1++;
  } elsif ($nextbest_d < $ratio*$best_d) {
    $n5++;
  } else {
    if ($nextbest_d > $x2) {
      $n2++;
    } elsif ($nextbest_d > $x1) {
      $n3++;
    } else {
      $n4++;
    }
  }
}
}
print "$n1, $n2, $n3, $n4, $n5 \n";
