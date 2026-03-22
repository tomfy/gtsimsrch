#!/usr/bin/perl -w
use strict;
use Scalar::Util qw(looks_like_number);

# usage:  quadrant_counts '3:0.04;16:0.3' <  <datafile>

# specify 2 cols with thresholds, e.g. '3:0.04;16:0.3'
# this defined 4 regions defined by
# whether the value is col 3 is > or < 0.04, and
# whether the values in col 16 is > or < 0.3
# e.g.
#    if col 3 <= 0.04, and col 16 < 0.3, then LL (for Lower Left) is appended to the end of the line
#    if col 3 > 0.04 and col 16 <= 0.3, the LR (for Lower Right) is appended to the end of the line
# If there are lines with one or both of these columns has a non-numerical value, it is not output.

my $qstring = shift // undef;
if (!defined $qstring) {
  print STDERR "# usage:  quadrant_counts '3:0.04;16:0.3' <  <datafile> \n";
  print STDERR "# must specify columns, thresholds and input file. Exiting \n";
  exit();
}
print "# $qstring \n";
print "# ", join(" ", @ARGV), "\n";
my $categorize = shift // 1;
my ($xstr, $ystr) = split(';', $qstring);
my ($xcol, $xthresh) = split(':', $xstr);
my ($ycol, $ythresh) = split(':', $ystr);
$xcol--;
$ycol--;
print STDERR "Xcol (zero-based): $xcol  Xthreshold: $xthresh\n";
print STDERR "Ycol (zero-based): $ycol  Ythreshold: $ythresh\n";

my ($LLcount, $ULcount, $LRcount, $URcount) = (0, 0, 0, 0);
my ($Lnan_count, $Unan_count, $nanL_count, $nanR_count) = (0, 0, 0, 0);
my ($nan_count, $nannan_count) = (0, 0);
while (<>) {
  next if(/^\s*#/);
  my @cols = split(" ", $_);
  my $xvalue = $cols[$xcol] // '-';
  my $yvalue = $cols[$ycol] // '-';
  # print stderr "$xvalue $yvalue   $xthresh $ythresh\n";
  #if ($xvalue =~ /\d+[.]?\d+/  and $yvalue =~ /\d+[.]?\d+/)
  #if(looks_like_number($xvalue) or looks_like_number($yvalue)){
  #$xvalue = '-' if($xcol < $#cols);
  #$yvalue = '-' if($ycol < $#cols);
  #print STDERR "$_" if(!looks_like_number($xvalue) or !looks_like_number($yvalue));
  chomp;
  if (! looks_like_number($xvalue)) {
    if (! looks_like_number($yvalue)) {
      $nannan_count++;
    } elsif ($yvalue <= $ythresh) {
      $Lnan_count++;
    } elsif ($yvalue > $ythresh) {
      $Unan_count++;
    }
  } elsif ($xvalue <= $xthresh) {
    if (! looks_like_number($yvalue)) {
      $nanL_count++;
    } elsif ($yvalue < $ythresh) {
      $LLcount++;
      print "$_ LL\n" if($categorize);
    } elsif ($yvalue > $ythresh) {
      $ULcount++;
      print "$_ UL\n" if($categorize);
    }
  } elsif ($xvalue > $xthresh) {      # x >= x threshold
    if (! looks_like_number($yvalue)) { #y nan
      $nanR_count++;
    } elsif ($yvalue <= $ythresh) {
      $LRcount++;
      print "$_ LR\n" if($categorize);
    } elsif ($yvalue > $ythresh) {
      $URcount++;
      print "$_ UR\n" if($categorize);
    }
  } 
}
my $left_count = $ULcount+$LLcount;
my $right_count = $URcount+$LRcount;
my $upper_count = $ULcount + $URcount;
my $lower_count = $LLcount + $LRcount;
my $total_2numbers_count = $upper_count + $lower_count;
$xcol++; $ycol++;		# back to unit-based col numbers.
printf("#\n");
print ("# Left  (Right) <-> col. $xcol  <=(>) $xthresh\n");
print ("# Lower (Upper) <-> col. $ycol  <=(>) $ythresh\n"); 
printf("#\n");
printf("#       |    Left  |  Right   |\n");
printf("# --------------------------------------\n");
printf("# Upper | %8d | %8d | %8d\n", $ULcount, $URcount, $upper_count);
printf("# --------------------------------------\n");
printf("# Lower | %8d | %8d | %8d\n", $LLcount, $LRcount, $lower_count);
printf("# --------------------------------------\n");
printf("#       | %8d | %8d | %8d\n", $left_count , $right_count, $total_2numbers_count);
printf("#\n");
printf("#(nan, Upper): %ld \n", $Unan_count);
printf("#(nan, Lower): %ld \n", $Lnan_count);
printf("#(Left, nan):  %ld \n", $nanL_count);
printf("#(Right, nan): %ld \n", $nanR_count);
printf("#(nan, nan):   %ld\n", $nannan_count);
