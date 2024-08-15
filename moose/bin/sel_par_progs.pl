#!/usr/bin/perl -w
use strict;

# read in pedigrees and hgmrs, d, z, etc, from pedigree test output file,
# and read in output from phased_parents.pl,
# and output lines with parent_pairs according to pedigree table


my %Fprog_parent_pairs = ();
my %Mprog_parent_pairs = ();
my $pedigree_test_output = shift;
open my $fhpt, "<", "$pedigree_test_output";
while (my $pt_line = <$fhpt>) {
  chomp $pt_line;
  my @cols = split(" ", $pt_line);
  my ($A, $F, $M) = @cols[0, 2, 3];
   $Fprog_parent_pairs{"$A $F"} .= "  $pt_line";
  $Mprog_parent_pairs{"$A $M"} .= "  $pt_line";
}
close $fhpt;

while (<>) {
  chomp;
  my @cols = split(" ", $_);
  my ($progeny, $parent) = @cols[0,1];
  print STDERR "$parent  $progeny \n";
  my $F_pt_line = $Fprog_parent_pairs{"$progeny $parent"} // "no F parent in pedigree ";
  print "$_ F  $F_pt_line\n";
  my $M_pt_line = $Mprog_parent_pairs{"$progeny $parent"} // "no M parent in pedigree ";
  print "$_ M  $M_pt_line\n";
}
