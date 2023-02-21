#!/usr/bin/perl -w
use strict;

my @marker_ids;
my @zeroes;
my @ones;
my @twos;
my @mds;
my $n_accessions = 0;
while(<>){
  next if(/^\s*#/);
  my @cols = split(" ", $_);
  if($cols[0] eq 'MARKER'){
    my $M = shift @cols;
    @marker_ids = @cols;
    my $N = scalar @marker_ids;
    @zeroes = (0) x $N;
    @ones = (0) x $N;
    @twos = (0) x $N;
    @mds = (0) x $N;
  }else{
    my $accid = shift @cols;
    while(my($i, $v) = each @cols){
      if($v eq 'NA'){
	$mds[$i]++;
      }elsif($v == 0){
	$zeroes[$i]++;
      }elsif($v == 1){
	$ones[$i]++;
      }elsif($v == 2){
	$twos[$i]++;
      }else{
	warn "dosage is $v ??\n";
      }
    }
  }
  $n_accessions++;
}

for(my $i=0; $i < scalar @marker_ids; $i++){
  my $minor_allele_count = #$twos[$i];
    ($zeroes[$i] < $twos[$i])? $zeroes[$i] : $twos[$i];
  my $dosage_sum = $ones[$i] + 2*$minor_allele_count;
  my $max_minor_alleles = 2*($n_accessions - $mds[$i]);
  my $maf_frequency = $dosage_sum/$max_minor_alleles;
  print $marker_ids[$i], "  ", $zeroes[$i], "  ", $ones[$i], "  ", $twos[$i], "  ", $mds[$i], "   $max_minor_alleles   $maf_frequency \n";
}
