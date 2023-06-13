#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw(min max sum);
use File::Spec 'splitpath';

my $reset_refs = 0;

GetOptions(
	   'Reset!' => \$reset_refs,
	   );

# test duplicatesearch comparing its output
# to known good output for several inputs and
# sets of run parameters.

my $failed_tests;
my $n_failed = 0;
# *************************************************************
#  yam941 data, unshuffled, est_max_dist = 0.2
my $test_name = "yam941 unshuffled est_max_dist=0.2";
my $input_file = "yam941_gp0.9.dsgs";
my $ref_stdout_file = "y941un.stdout";
my $test_stdout_file = "y941un.stdout_x";
my $ref_dsout_file = "y941_gp0.9_maf0.08_maxd0.2_md0.25_un.dsout";
my $test_dsout_file = "y941_gp0.9_maf0.08_maxd0.2_md0.25_un.dsout_x";

if($reset_refs){
  reset_refs($input_file, $ref_stdout_file, $ref_dsout_file, " -est_max_dist 0.2  -un ");
}else{
  $failed_tests = test($test_name, $input_file, $ref_stdout_file, $test_stdout_file, $ref_dsout_file, $test_dsout_file, " -est_max_dist 0.2  -un ");
  $n_failed++ if($failed_tests);
}


#**************************************************************
#  yam941 data, shuffled, est_max_dist = 0.2
$test_name = "yam941 shuffled est_max_dist=0.2";
$ref_stdout_file = "y941.stdout";
$test_stdout_file = "y941.stdout_x";

$ref_dsout_file = "y941_gp0.9_maf0.08_maxd0.2_md0.25.dsout";
$test_dsout_file = "y941_gp0.9_maf0.08_maxd0.2_md0.25.dsout_x";

if($reset_refs){
  reset_refs($input_file, $ref_stdout_file, $ref_dsout_file, " -est_max_dist 0.2 ");
}else{
  $failed_tests = test($test_name, $input_file, $ref_stdout_file, $test_stdout_file, $ref_dsout_file, $test_dsout_file, " -est_max_dist 0.2 ");
    $n_failed++ if($failed_tests);

  #print "Passed ", $n_tests-$failed_test, " of $n_tests.\n";
}

#**************************************************************
# cassava, shuffled, est_max_dist = 0.17

$test_name = "cassava shuffled est_max_dist=0.17";
$input_file = "cassava_gp0.9.dsgs";

$ref_stdout_file = "cassava.stdout";
$test_stdout_file = "cassava.stdout_x";

$ref_dsout_file = "cassava_gp0.9_maf0.08_maxd0.17_md0.25.dsout";
$test_dsout_file = "cassava_gp0.9_maf0.08_maxd0.17_md0.25.dsout_x";

if($reset_refs){
  reset_refs($input_file, $ref_stdout_file, $ref_dsout_file, " -est_max_dist 0.17 ");
}else{
  $failed_tests = test($test_name, $input_file, $ref_stdout_file, $test_stdout_file, $ref_dsout_file, $test_dsout_file, " -est_max_dist 0.17 ");
    $n_failed++ if($failed_tests);

  #print "Passed ", $n_tests-$failed_test, " of $n_tests.\n";
}

#**************************************************************
# cassava, shuffled, est_max_dist based on random sample of distances

$test_name = "cassava_shuffled est_max_dist=auto";
$input_file = "cassava_gp0.9.dsgs";

$ref_stdout_file = "cassava_a.stdout";
$test_stdout_file = "cassava_a.stdout_x";

$ref_dsout_file = "cassava_gp0.9_maf0.08_md0.25.dsout";
$test_dsout_file = "cassava_gp0.9_maf0.08_md0.25.dsout_x";

if($reset_refs){
  reset_refs($input_file, $ref_stdout_file, $ref_dsout_file, "");
}else{
  $failed_tests = test($test_name, $input_file, $ref_stdout_file, $test_stdout_file, $ref_dsout_file, $test_dsout_file, "");
  $n_failed++ if($failed_tests);

  #print "Passed ", $n_tests-$failed_test, " of $n_tests.\n";
}

##################################################################################

sub reset_refs{
  my $input_file = shift;
  my $ref_stdout_file = shift;
  my $ref_dsout_file = shift;
  my $extra_arguments = shift;

  my $stderr_file = "stderr";

  system "chmod 666 $ref_stdout_file";
  system "chmod 666 $ref_dsout_file";

  my $command = "duplicatesearch  -in $input_file  -out $ref_dsout_file  2> $stderr_file  > $ref_stdout_file ";
  $command .= " -maf 0.08  -accession_max 0.5  -seed 13572467  -t 4 ";
  $command .= $extra_arguments;

  print "$command \n";
  system "$command";

  system "chmod 444 $ref_stdout_file";
  system "chmod 444 $ref_dsout_file";

  print "$ref_stdout_file and $ref_dsout_file reset.\n";
}

sub test{
  my $test_name = shift;
  my $input_file = shift;
  my $ref_stdout_file = shift;
  my $test_stdout_file = shift;
  my $ref_dsout_file = shift;
  my $test_dsout_file = shift;
  my $extra_arguments = shift; # e.g. " -un ";

  my $stderr_file = "stderr";
  my $failed_tests = "";

  my $command = "duplicatesearch  -in $input_file  -out $test_dsout_file  2> $stderr_file  > $test_stdout_file ";
  $command .= " -maf 0.08  -accession_max 0.5  -seed 13572467  -t 4 ";
  $command .= $extra_arguments;

  system "$command";

  my $bad_stderr = 0;
  my $check_genotypes_grep  = `grep 'Successfully completed check_genotypesset' $stderr_file`;
  $bad_stderr = 1 if(length $check_genotypes_grep == 0);
  my $stderr_wc_str = `wc -l $stderr_file`;
  $bad_stderr = 1 if(! ($stderr_wc_str =~ /^\s*6\s/));
  $failed_tests .= "stderr " if($bad_stderr > 0);

  my $diffstr =  `diff $ref_stdout_file $test_stdout_file`;
  my $nnl = () = $diffstr =~ /\n/g;
  #printf("Past test: diff wrt $ref_stdout_file. diff has %d lines.\n", $nnl);
  $failed_tests .=  "stdout " if($nnl != 4);
  # print $diffstr;

  $diffstr = `diff $ref_dsout_file $test_dsout_file`;
  $nnl = () = $diffstr =~ /\n/g;
  #printf("Past test: diff wrt $ref_dsout_file.  diff has %d lines.\n", $nnl);
  $failed_tests .= "dsout " if($nnl != 4);
  # print "[" . $diffstr . "]" . "\n" if($nnl > 0);
  unlink ($test_stdout_file, $test_dsout_file) if($failed_tests);
  print "$test_name : ", ($failed_tests)? "failed\n" : "passed\n";
  return $failed_tests;
}
