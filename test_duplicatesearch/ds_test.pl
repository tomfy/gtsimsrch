#!/usr/bin/perl -w
use strict;

my $reset_refs = shift // 0;

# test duplicatesearch comparing its output
# to known good output for several inputs and
# sets of run parameters.

my ($n_tests, $n_failed_tests);
# *************************************************************
#  yam941 data, unshuffled, max_dist = 0.2

my $input_file = "yam941_gp0.9.dsgs";
my $ref_stdout_file = "y941un.stdout";
my $test_stdout_file = "y941un.stdout_x";
my $ref_dsout_file = "y941_gp0.9_maf0.08_maxd0.2_md0.25_un.dsout";
my $test_dsout_file = "y941_gp0.9_maf0.08_maxd0.2_md0.25_un.dsout_x";

if($reset_refs){
  resetx($input_file, $ref_stdout_file, $ref_dsout_file, " -max_dist 0.2  -un ");
}else{
  ($n_tests, $n_failed_tests) = test($input_file, $ref_stdout_file, $test_stdout_file, $ref_dsout_file, $test_dsout_file, " -est_max_dist 0.2  -un ");
  print "Passed ", $n_tests-$n_failed_tests, " of $n_tests.\n";
}

#**************************************************************
#  yam941 data, shuffled, max_dist = 0.2

$ref_stdout_file = "y941.stdout";
$test_stdout_file = "y941.stdout_x";

$ref_dsout_file = "y941_gp0.9_maf0.08_maxd0.2_md0.25.dsout";
$test_dsout_file = "y941_gp0.9_maf0.08_maxd0.2_md0.25.dsout_x";

if($reset_refs){
  resetx($input_file, $ref_stdout_file, $ref_dsout_file, " -max_dist 0.2 ");
}else{
  ($n_tests, $n_failed_tests) = test($input_file, $ref_stdout_file, $test_stdout_file, $ref_dsout_file, $test_dsout_file, " -est_max_dist 0.2 ");
  print "Passed ", $n_tests-$n_failed_tests, " of $n_tests.\n";
}

#**************************************************************
# cassava, shuffled, max_dist = 0.17

$input_file = "cassava_gp0.9.dsgs";

$ref_stdout_file = "cassava.stdout";
$test_stdout_file = "cassava.stdout_x";

$ref_dsout_file = "cassava_gp0.9_maf0.08_maxd0.17_md0.25.dsout";
$test_dsout_file = "cassava_gp0.9_maf0.08_maxd0.17_md0.25.dsout_x";

if($reset_refs){
  resetx($input_file, $ref_stdout_file, $ref_dsout_file, " -max_dist 0.17 ");
}else{
  ($n_tests, $n_failed_tests) = test($input_file, $ref_stdout_file, $test_stdout_file, $ref_dsout_file, $test_dsout_file, " -est_max_dist 0.17 ");
  print "Passed ", $n_tests-$n_failed_tests, " of $n_tests.\n";
}

#**************************************************************
# cassava, shuffled, max_dist based on random sample of distances

$input_file = "cassava_gp0.9.dsgs";

$ref_stdout_file = "cassava_a.stdout";
$test_stdout_file = "cassava_a.stdout_x";

$ref_dsout_file = "cassava_gp0.9_maf0.08_md0.25.dsout";
$test_dsout_file = "cassava_gp0.9_maf0.08_md0.25.dsout_x";

if($reset_refs){
  resetx($input_file, $ref_stdout_file, $ref_dsout_file, "");
}else{
  ($n_tests, $n_failed_tests) = test($input_file, $ref_stdout_file, $test_stdout_file, $ref_dsout_file, $test_dsout_file, "");
  print "Passed ", $n_tests-$n_failed_tests, " of $n_tests.\n";
}

sub resetx{
  my $input_file = shift;
  my $ref_stdout_file = shift;
  my $ref_dsout_file = shift;
  my $extra_arguments = shift;

  my $stderr_file = "stderr";

  my $command = "duplicatesearch  -in $input_file  -out $ref_dsout_file  2> $stderr_file  > $ref_stdout_file ";
  $command .= " -maf 0.08  -accession_max 0.5  -seed 13572467  -t 4 ";
  $command .= $extra_arguments;

  system "$command";
  print "$ref_stdout_file and $ref_dsout_file reset.\n";
}

sub test{
  my $input_file = shift;
  my $ref_stdout_file = shift;
  my $test_stdout_file = shift;
  my $ref_dsout_file = shift;
  my $test_dsout_file = shift;
  my $extra_arguments = shift; # e.g. " -un ";

  my $stderr_file = "stderr";
  my $n_failed_tests = 0;

  my $command = "duplicatesearch  -in $input_file  -out $test_dsout_file  2> $stderr_file  > $test_stdout_file ";
  $command .= " -maf 0.08  -accession_max 0.5  -seed 13572467  -t 4 ";
  $command .= $extra_arguments;
  
  #print "$command \n";
  #print `ls -lt *.dsout*`, "\n";

  system "$command";

  my $diffstr =  `diff $ref_stdout_file $test_stdout_file`;
  my $nnl = () = $diffstr =~ /\n/g;
  printf("Past test: diff wrt $ref_stdout_file. diff has %d lines.\n", $nnl);
  $n_failed_tests++ if($nnl != 4);
  # print $diffstr;

  $diffstr = `diff $ref_dsout_file $test_dsout_file`;
  $nnl = () = $diffstr =~ /\n/g;
  printf("Past test: diff wrt $ref_dsout_file.  diff has %d lines.\n", $nnl);
  $n_failed_tests++ if($nnl != 4);
  # print "[" . $diffstr . "]" . "\n" if($nnl > 0);
  unlink ($test_stdout_file, $test_dsout_file, $stderr_file) if($n_failed_tests == 0);
  return (2, $n_failed_tests);
}
