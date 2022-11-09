#!/usr/bin/perl -w
use strict;

my $the_col = shift // 18;	# unit-based
$the_col--;			# now zero-based

my %accid_solns = ();
while (my $line = <>) {
  next if($line =~ /^\s*#/);
   my @cols = split(" ", $line);
   my $value = $cols[$the_col];
  my $accid = shift @cols;
  my $parents = $cols[0] . " " . $cols[1];
  my $cross_type = ($cols[0] eq $cols[1])? "self" : "bip";
  my $x = join(" ", @cols);
  my $vc = $value . "\t" .  "$cross_type $x ";
  if (exists $accid_solns{$accid}) {
    push @{$accid_solns{$accid}}, $vc; # "$parents  $value"; 
  } else {
    $accid_solns{$accid} = [$vc];
  }
}

for my $anid (keys %accid_solns) {
  my $solutions = $accid_solns{$anid};
  my @sorted_solns = sort {val($a) <=> val($b)} @$solutions;
  print "$anid  ", $sorted_solns[0];
  print "  ", $sorted_solns[1] if(scalar @sorted_solns > 1);
  print "\n"; # join("  ", @sorted_solns), "\n";
}


sub val{
  my $vc = shift;
  my ($value, $xx) = split("\t", $vc);
  return $value;
}
