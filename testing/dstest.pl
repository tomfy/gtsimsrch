#!/usr/bin/perl -w
use strict;
# use File::Spec qw(splitpath);
use Test::More tests => 1;
use File::Basename 'dirname';
use List::Util qw 'min max sum';

use Cwd 'abs_path';
my ( $bindir, $libdir, $datadir );
BEGIN {     # this has to go in Begin block so happens at compile time
  $bindir =
    dirname( abs_path(__FILE__) ) ; # the directory containing this script
  $libdir = $bindir . '/../moose/lib';
  $libdir = abs_path($libdir);	# collapses the bin/../lib to just lib
  $datadir = $bindir; # . '/../test_data';
  $datadir = abs_path($datadir);	# collapses the bin/../data to just data
}
use lib $libdir;

my $dosages_file = $datadir . '/8698x5621.dsgs';
my $std_dsout_file = $datadir . '/8698x5621.dsout_std';
# print "# dosages file: $dosages_file\n";
# print "# duplicate_search output standard file: $std_dsout_file \n";
my $test_dsout_file = '8698x5621.dsout_test';

my $rng_seed = 123456789;

my $ds_command = "duplicate_search  -in $dosages_file  -out $test_dsout_file  -seed $rng_seed  > /dev/null  2> /dev/null";
# print "$ds_command \n";
system "$ds_command";

# print "diff $std_dsout_file $test_dsout_file \n";
my $diff_out = `diff $std_dsout_file $test_dsout_file`;

print "diff output:[$diff_out]\n";
ok($diff_out eq '', 'diff_of_outputs');

unlink $test_dsout_file;
