#!/usr/bin/perl -w
use strict;

# input is a vcf file with GT and AD fields
my $ploidy = shift // 6;

while (my $line = <>) {
  next if($line =~ /^\s*#/);
  my @cols = split(" ", $line);
  my $fstr = $cols[8];
  die if(! $fstr =~ /^GT\:AD\:DP$/);
  @cols = @cols[9..$#cols];

  for my $c (@cols) {
    my ($gt, $ad, $dp) = split(":", $c);
    #  print "$gt    $ad    $dp    ";
    if ($gt =~ /\./) {
      print "x   x  ";
    } else {
      my $ref_count = $gt =~ tr/0/0/;
      my $alt_count = $gt =~ tr/1/1/;
      print "$alt_count  $dp  ";
    }
    my ($ref_depth, $alt_depth) = split(",", $ad);
    die if($ref_depth + $alt_depth != $dp);
    my $s = "$ref_depth  $alt_depth   ";
    $s .= ($dp > 0)? $ploidy*$alt_depth/$dp . "\n" : 2*$ploidy . "\n";
    print $s;
  }
}
    
  
