#!/usr/bin/perl -w
use strict;

# test duplicatesearch comparing its output
# to known good output for several inputs and
# sets of run parameters.

my $input_file = "yam941_gp0.9.dsgs";

my $ref_stdout_file = "y941un.stdout";
my $test_stdout_file = "y941un.stdout_x";

my $ref_dsout_file = "y941_gp0.9_maf0.08_maxd0.2_md0.25_un.dsout";
my $test_dsout_file = "y941_gp0.9_maf0.08_maxd0.2_md0.25_un.dsout_x";

my $command = "duplicatesearch -in $input_file  -out $test_dsout_file  2> stderr  > $test_stdout_file ";
  $command .= " -maf 0.08  -max_dist 0.2  -accession_max 0.5  -seed 13572467  -marker_max 0.25 ";
  $command .= " -un  -t 4";

system "$command";

system "diff $ref_stdout_file $test_stdout_file";
print "Past test 1a: diff wrt $ref_stdout_file\n";

system "diff $ref_dsout_file $test_dsout_file";
print "Past test 1b: diff wrt $ref_dsout_file\n";

#**************************************************************
$ref_stdout_file = "y941.stdout";
$test_stdout_file = "y941.stdout_x";

$ref_dsout_file = "y941_gp0.9_maf0.08_maxd0.2_md0.25.dsout";
$test_dsout_file = "y941_gp0.9_maf0.08_maxd0.2_md0.25.dsout_x";

 $command = "duplicatesearch -in $input_file  -out $test_dsout_file  2> stderr  > $test_stdout_file ";
  $command .= " -maf 0.08  -max_dist 0.2  -accession_max 0.5  -seed 13572467  -marker_max 0.25 ";
  $command .= " -t 4";


system "$command";

system "diff $ref_stdout_file $test_stdout_file";
print "Past test 2a: diff wrt $ref_stdout_file\n";

system "diff $ref_dsout_file $test_dsout_file";
print "Past test 2b: diff wrt $ref_dsout_file\n";
