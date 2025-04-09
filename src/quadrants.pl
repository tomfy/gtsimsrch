#!/usr/bin/perl -w
use strict;

# usage:  quadrant_counts '3:0.04;16:0.3' <  <datafile>

# specify 2 cols with thresholds, e.g. '3:0.04;16:0.3'
# this defined 4 regions defined by
# whether the value is col 3 is > or < 0.04, and
# whether the values in col 16 is > or < 0.3


my $qstring = shift // undef;
if (!defined $qstring) {
  print STDERR "# usage:  quadrant_counts '3:0.04;16:0.3' <  <datafile> \n";
  print STDERR "# must specify columns, thresholds and input file. Exiting \n";
  exit();
}
my $categorize = shift // 1;
my ($xstr, $ystr) = split(';', $qstring);
my ($xcol, $xthresh) = split(':', $xstr);
my ($ycol, $ythresh) = split(':', $ystr);
$xcol--;
$ycol--;
print STDERR "Xcol (zero-based): $xcol  Xthreshold: $xthresh\n";
print STDERR "Ycol (zero-based): $ycol  Ythreshold: $ythresh\n";

my ($LLcount, $ULcount, $LRcount, $URcount) = (0, 0, 0, 0);
my $nan_count = 0;
while (<>) {
  next if(/^\s*#/);
  my @cols = split(" ", $_);
  my $xvalue = $cols[$xcol];
  my $yvalue = $cols[$ycol];
  if ($xvalue =~ /\d+[.]\d+/  and $yvalue =~ /\d+[.]\d+/) {
    chomp;
    if ($xvalue < $xthresh) {
      if ($yvalue < $ythresh) {
	$LLcount++;
	print "$_ LL\n" if($categorize);
      } else {
	$ULcount++;
	print "$_ UL\n" if($categorize);
      }
    } else {			# x >= x threshold
      if ($yvalue < $ythresh) {
	$LRcount++;
	print "$_ LR\n" if($categorize);
      } else {
	$URcount++;
	print "$_ UR\n" if($categorize);
      }
    }
  } else {
    $nan_count++;
  }
}
printf("# count of accessions with <2 values: %8d\n", $nan_count);
printf("# %8d %8d\n", $ULcount, $URcount);
printf("# %8d %8d\n", $LLcount, $LRcount);
