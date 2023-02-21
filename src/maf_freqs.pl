#!/usr/bin/perl -w
use strict;

my $rectify = shift // 1;
my $acc_max_md_fraction = shift // 0.5;
my $marker_max_md_fraction = shift // 0.25;
my @marker_ids;
my @zeroes;
my @ones;
my @twos;
my @mds;
my $N;
my $n_accessions = 0;
my $n_ok_accessions = 0;
while (<>) {
  next if(/^\s*#/);
  my @cols;
  if ( /^\s*MARKER/ ) {		#$cols[0] eq 'MARKER'){
    @cols = split(" ", $_);
    my $M = shift @cols;
    @marker_ids = @cols;
    $N = scalar @marker_ids;
    @zeroes = (0) x $N;
    @ones = (0) x $N;
    @twos = (0) x $N;
    @mds = (0) x $N;
  } else {
    s/^\s*(\S+)//;
    my $accid = $1;
    my $NA_count = () = /NA/gi;
    $n_accessions++;
    next if($NA_count > $acc_max_md_fraction*$N); # accession has excessive missing data
    @cols = split(" ", $_);
      while (my($i, $v) = each @cols) {
	if ($v eq 'NA') {
	  $mds[$i]++;
	} elsif ($v == 0) {
	  $zeroes[$i]++;
	} elsif ($v == 1) {
	  $ones[$i]++;
	} elsif ($v == 2) {
	  $twos[$i]++;
	} else {
	  warn "dosage is $v ??\n";
	}
      }
  }
  $n_ok_accessions++;
}

for (my $i=0; $i < scalar @marker_ids; $i++) {
  my $minor_allele_count = ($rectify)?  (($zeroes[$i] < $twos[$i])? $zeroes[$i] : $twos[$i]) : $twos[$i];
   
  my $dosage_sum = $ones[$i] + 2*$minor_allele_count;
  next if($dosage_sum == 0);
  next if($mds[$i] > $marker_max_md_fraction*$n_ok_accessions);
  my $total_observed_alleles = 2*($n_ok_accessions - $mds[$i]);
  my $maf_frequency = $dosage_sum/$total_observed_alleles;
  
  print $marker_ids[$i], "  ", $zeroes[$i], "  ", $ones[$i], "  ", $twos[$i], "  ", $mds[$i], "   $total_observed_alleles   $maf_frequency \n";
}
print "# n_accessions: $n_accessions  n_ok_accessions: $n_ok_accessions\n";
