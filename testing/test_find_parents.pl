#!/usr/bin/perl -w
use strict;

use File::Basename 'dirname';
use List::Util qw 'min max sum';
use Capture::Tiny qw 'capture'; 

use Cwd 'abs_path';
my ( $bindir, $libdir, $datadir );
BEGIN {     # this has to go in Begin block so happens at compile time
  $bindir =
    dirname( abs_path(__FILE__) ) ; # the directory containing this script
  $libdir = $bindir . '/../moose/lib';
  $libdir = abs_path($libdir);	# collapses the bin/../lib to just lib
  $datadir = $bindir;		# . '/../test_data';
  $datadir = abs_path($datadir); # collapses the bin/../data to just data
}
use lib $libdir;

my $tests_to_run = shift // '0'; # e.g.: '0,2,3,4' or all
if ($tests_to_run eq 'all') {
  $tests_to_run = '0,2,3,4';
}
my @tsts = split(",", $tests_to_run);
my %tests2run = ();
for (@tsts) {
  $tests2run{$_} = 1;
}
my $n_tests_to_run = scalar keys %tests2run;
use Test::More;
# test find_parents in 5 modes:
# no pedigree file;
# check pedigrees from file and
#   0) don't search for alternative pedigrees;
#   1) search, consider as parents only accessions which are parents in pedigree file;
#   2) search if pedigree in file is bad, or nonexistent;
#   3) search for alternatives for all accessions.

if (exists $tests2run{0}) {
  # just check pedigrees in file (i.e. -alt 0 (which is default))
  my $output_standard_file = $datadir . "/alt0.fpout_std";
  my $output_test_file = "alt0.fpout";
  my $dosage_input_file = $datadir . "/u_wchroms.dsgs";
  my $pedigree_input_file = $datadir . "/u_ped_table";
  my $fp_command = "find_parents -seed 1234567 -in $dosage_input_file -ped $pedigree_input_file -out $output_test_file";
  my ($stdout, $stderr) = capture{ system "$fp_command" };
  my $diff_out = `diff $output_standard_file $output_test_file`;
  #unlink $output_test_file;

  ok($diff_out eq '', 'alt0  output agrees with standard.');
}

if (exists $tests2run{1}) {
  # check pedigrees in file and find alternatives,
  # but only consider as possible parents those which are parents (of something) according to pedigree file.
  # (i.e. -alt 1 (which is default))
  my $output_standard_file = $datadir . "/alt1.fpout_std";
  my $output_test_file = "alt1.fpout";
  my $dosage_input_file = $datadir . "/u_wchroms.dsgs";
  my $pedigree_input_file = $datadir . "/u_ped_table";
  my $fp_command = "find_parents -seed 1234567 -in $dosage_input_file -ped $pedigree_input_file -alt 1 -out $output_test_file";
  my ($stdout, $stderr) = capture{ system "$fp_command" };
  my $diff_out = `diff $output_standard_file $output_test_file`;
  # unlink $output_test_file;

  ok($diff_out eq '', 'alt1  output agrees with standard.');
}

# check pedigrees in file and, if bad, search for alternatives/
# (i.e. -alt 2)
if (exists $tests2run{2}) {
  my $output_standard_file = $datadir . "/alt2.fpout_std";
  my $output_test_file = "alt2.fpout";
  my $dosage_input_file = $datadir . "/u_wchroms.dsgs";
  my $pedigree_input_file = $datadir . "/u_ped_table";
  my $fp_command = "find_parents -seed 1234567 -in $dosage_input_file -ped $pedigree_input_file -alt 2 -out $output_test_file";
  my ($stdout, $stderr) = capture{ system "$fp_command" };
  my $diff_out = `diff $output_standard_file $output_test_file`;
  unlink $output_test_file;

  ok($diff_out eq '', 'alt2  output agrees with standard.'); 
}

# check pedigrees in file and, if bad, search for alternatives/
#(i.e. -alt 3)
if (exists $tests2run{3}) {
  my $output_standard_file = $datadir . "/alt3.fpout_std";
  my $output_test_file = "alt3.fpout";
  my $dosage_input_file = $datadir . "/u_wchroms.dsgs";
  my $pedigree_input_file = $datadir . "/u_ped_table";
  my $fp_command = "find_parents -seed 1234567 -in $dosage_input_file -ped $pedigree_input_file -alt 3 -out $output_test_file";
  my ($stdout, $stderr) = capture{ system "$fp_command"};
  my $diff_out = `diff $output_standard_file $output_test_file`;
  unlink $output_test_file; 

  ok($diff_out eq '', 'alt3  output agrees with standard.');
}

# no-pedigree case:
if (exists $tests2run{4}) {
  my $output_standard_file = $datadir . "/noped.fpout_std";
  my $output_test_file = 'noped.fpout_test';
  my $dosage_input_file = $datadir . "/u_wchroms.dsgs";
  my $fp_command = "find_parents -seed 1234567 -in $dosage_input_file -out $output_test_file";
  my ($stdout, $stderr) = capture{ system "$fp_command" };
  my $diff_out = `diff $output_standard_file $output_test_file`;
  unlink $output_test_file;
  
  ok($diff_out eq '', 'no-pedigree case agrees with standard.'); 
}
done_testing( $n_tests_to_run );

