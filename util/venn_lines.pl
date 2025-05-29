#!/usr/bin/perl -w
use strict;
use List::Util qw(min max sum);

# useful for seeing whether only difference between files in that lines
# are in a different order.
# default is to skip comments (lines starting with # )
# outputs count of lines in file 1 only, in both files, and in file 2 only.

# usage example:
# vennlines.pl file1 file2      # will ignore comment lines
# vennlines.pl file1 file2  0    # also compare comment lines
# vennlines.pl file1 file2  1 1    # skip comments, but output lines found in only one of the files
{
  my %file1lines = ();
  my %file2lines = ();

  my $file1 = shift;
  my $file2 = shift;
  my $skip_comments = shift // 1;
  my $verbose = shift // 0;

  store($file1, \%file1lines, $skip_comments);
  store($file2, \%file2lines, $skip_comments);

  my @one_onlies = ();
  my @boths = ();
  my @two_onlies = ();

  while (my ($line, $n) = each %file1lines) {
    if (exists $file2lines{$line}) {
      push @boths, $line;
    } else {
      push @one_onlies, $line;
    }
  }

  my $both_count = 0;
  while (my ($line, $n) = each %file2lines) {
    if (exists $file1lines{$line}) {
      $both_count++;
    } else {
      push @two_onlies, $line;
    }
  }
  #print scalar keys %file1lines, "   ", scalar keys %file2lines, "   ", sum(values %file1lines), "  ", sum(values %file2lines), "\n";
  print "count of lines in file 1 only: ", scalar @one_onlies, "\n";
  print "count of lines in both files:  ", $both_count, "\n";
  print "count of lines in both files:  ", scalar @boths, " (as a check, should be same as previous line)\n";
  print "count of lines in file 2 only: ", scalar @two_onlies, "\n";
  if($verbose){
    print "lines in $file1 only: \n";
    print join("", @one_onlies);
    print "lines in $file2 only: \n";
    print join("", @two_onlies);
  }
}

sub store{
  my $file = shift;
  my $filelines = shift;	# hashref
  my $skip_comments = shift;

  open my $fh, "<", $file or die "couldn't open $file for reading.\n";
  while (<$fh>) {
    next if($skip_comments  and  /^\s*#/); # skip comments
    if (/^\s*(\S+)/) {
      $filelines->{$_}++;
    }
  }
  return $filelines;
  close $fh;
}
