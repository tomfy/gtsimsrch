#!/usr/bin/perl -w
use strict;

my $the_col = shift // 16;	# unit-based
$the_col--;			# now zero-based

my %accid_solns = ();
while (my $line = <>) {
  next if($line =~ /^\s*#/);
  my @cols = split(" ", $line);
  my $accid = $cols[0];
  my $parents = $cols[1] . " " . $cols[2];
  my $value = $cols[$the_col];
  if (exists $accid_solns{$accid}) {
    push @{$accid_solns{$accid}}, $value; # "$parents  $value"; 
  } else {
    $accid_solns{$accid} = [$value];
  }
}

for my $anid (keys %accid_solns) {
  my $solutions = $accid_solns{$anid};
  my @sorted_solns = sort {$a <=> $b} @$solutions;
  print "$anid  ", join("  ", @sorted_solns), "\n";
}
