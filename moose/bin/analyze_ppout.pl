#!/usr/bin/perl -w
use strict;
use Getopt::Long;

# analyze output of phased_parents.pl
# each input line has information about 1 accession pair and 1 pair of homologous chromosomes.
# for each, choose the chromosome of the homologous pair in the 'offspring' which minimizes required crossovers.
# and sum over the chromosome, also counting the number of heterozygous gts in the parent.

my $pp_out_file = undef; # which is input to analyze_ppout
my $output_file = "app.out";

GetOptions(
	   'phased_parents_output_file|pp_out_file=s' => \$pp_out_file,
	   'output_file=s' => \$output_file,
	  );

if(!defined $pp_out_file){
  print_help();
  exit;
}

my %pp_info = (); # keys: parent-progeny id pairs (e.g.  "parid progid"

open my $fhin, "<", "$pp_out_file";
while(my $pp_line = <$fhin>){
  next if($pp_line =~ /^\s*(#|$)/);
  my ($progid, $parid, $i_chrom, $hgmr, $XA, $XB, $Xmin, $Xmax,
      # $ndubA, $ndubB,
      $parent_het_count, $type, $direction) = split(" ", $pp_line);
  my $ppid = "$progid $parid";
  if(!exists $pp_info{$ppid}){
    my $info = {'xmin' => $Xmin, 'xmax' => $Xmax, 'parent_hetzygs' => $parent_het_count, 'hgmr' => $hgmr, 'type' => $type}; 
    $pp_info{$ppid} = $info;
  }else{
    $pp_info{$ppid}->{xmin} += $Xmin;
    $pp_info{$ppid}->{xmax} += $Xmax;
    $pp_info{$ppid}->{parent_hetzygs} += $parent_het_count;
  }
}
close $fhin;

open my $fhout, ">", "$output_file";
while (my($pp, $info) = each %pp_info) {
  my ($xmin, $xmax, $parhets, $hgmr, $type) = map($info->{$_}, ('xmin', 'xmax', 'parent_hetzygs', 'hgmr', 'type'));
  print STDERR "$pp   $xmin  $xmax  $parhets  $hgmr\n";
  my $xmin_norm = ($parhets > 0)? $xmin/$parhets : '-';
  my $xmax_norm = ($parhets > 0)? $xmax/$parhets : '-';

  my ($p1, $p2) = split(" ", $pp);
  my $rev_pp = "$p2 $p1";
  my $rev_info = $pp_info{$rev_pp} // undef;
   print $fhout "$pp  $type $hgmr  ";
  if (defined $rev_info) {
    my ($rev_xmin, $rev_xmax, $rev_parhets, $rev_hgmr) = map($rev_info->{$_}, ('xmin', 'xmax', 'parent_hetzygs', 'hgmr'));
    my $rev_xmin_norm = $rev_xmin/$rev_parhets // '-';
    my $rev_xmax_norm = $rev_xmax/$rev_parhets // '-';
      print $fhout "$xmin $xmax  $parhets  $xmin_norm $xmax_norm  ",
      "$rev_xmin $rev_xmax  $rev_parhets  $rev_xmin_norm $rev_xmax_norm \n";
  } else {
      print $fhout "$xmin $xmax  $parhets  $xmin_norm $xmax_norm \n";
  }
}
close $fhout;

sub print_help{
  print STDOUT "usage: analyze_ppout -pp_out_file  <output file from phased_parents> [-out <output filename>]\n";
}
