#!/usr/bin/perl -w
use strict;

# look at pedigree_test output,
# there may be, due to uniquification (or otherwise?)
# accessions for which more than one pedigree is given (on separate lines)
# 

my %progid_lines = ();
while (<>) {
  next if(/^\s*#/);
  my @cols = split(" ", $_);
  my $progid = $cols[0];
  my $parpair = $cols[2] . " " . $cols[3];
  my $D = $cols[17];
  my $xhgmrF = $cols[19];
  my $xhgmrM = $cols[20];
  my $X = "$parpair  $D";
  if (!exists $progid_lines{$progid}) { # first line with this progeny id
    $progid_lines{$progid} = [$X];
  } else {
    push @{$progid_lines{$progid}}, $X;
  }
}




while (my($id, $count) = each %progid_lines) {
  my @the_Xs = @{$progid_lines{$id}};
  if (scalar @the_Xs > 1) {
    @the_Xs = sort {DfromX($a) <=> DfromX($b)} @the_Xs;
    print STDERR "$id  ", join("  ", @the_Xs), "\n";
  }
}

sub DfromX{
  my $X = shift;
  return ($X =~ /(\S+)\s*$/)? $1 : undef;
}

# sub S{
#   my $s1 = shift;
#   my $s2 = shift;
#   my @cols1 = split(" ", $s1);
#   my @cols2 = split(" ", $s2);
#   return $cols1[-1] <=> $cols2[-1];
# }
