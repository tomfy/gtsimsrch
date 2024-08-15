#!/usr/bin/perl -w
use strict;

# analyze output of phased_parents.pl
# each input line has information about 1 accession pair and 1 pair of homologous chromosomes.
# for each, choose the chromosome of the pair which minimized required recombinations.

my %pp_info = (); # keys: parent-progeny id pairs (e.g.  "parid progid"

while(<>){
  next if(/^\s*#/);
	my ($i, $j, $parid, $progid, $i_chrom, $hgmr, $XA, $XB, $Xmin, $Xmax, $parent_het_count, $pFM) = split(" ", $_);
  my $ppid = "$parid $progid";
  if(!exists $pp_info{$ppid}){
    my $info = {'xmin' => $Xmin, 'xmax' => $Xmax, 'parent_hetzygs' => $parent_het_count, 'hgmr' => $hgmr}; 
    $pp_info{$ppid} = $info;
  }else{
    $pp_info{$ppid}->{xmin} += $Xmin;
    $pp_info{$ppid}->{xmax} += $Xmax;
    $pp_info{$ppid}->{parent_hetzygs} += $parent_het_count;
  }
}

while (my($pp, $info) = each %pp_info) {
  my ($xmin, $xmax, $parhets, $hgmr) = map($info->{$_}, ('xmin', 'xmax', 'parent_hetzygs', 'hgmr'));
  print STDERR "$pp   $xmin  $xmax  $parhets  $hgmr\n";
  my $xmin_norm = ($parhets > 0)? $xmin/$parhets : -1;
  my $xmax_norm = ($parhets > 0)? $xmax/$parhets : -1;

  my ($p1, $p2) = split(" ", $pp);
  my $rev_pp = "$p2 $p1";
  my $rev_info = $pp_info{$rev_pp} // undef;
  if (defined $rev_info) {
    my ($rev_xmin, $rev_xmax, $rev_parhets, $rev_hgmr) = map($rev_info->{$_}, ('xmin', 'xmax', 'parent_hetzygs', 'hgmr'));
    my $rev_xmin_norm = $rev_xmin/$rev_parhets // -1;
    my $rev_xmax_norm = $rev_xmax/$rev_parhets // -1;
    print "$pp  $hgmr  ",
      "$xmin $xmax  $parhets  $xmin_norm $xmax_norm  ",
      "$rev_xmin $rev_xmax  $rev_parhets  $rev_xmin_norm $rev_xmax_norm ",
      "\n";
  } else {
    print "$pp  $hgmr  ",
      "$xmin $xmax  $parhets  $xmin_norm $xmax_norm  ",
      # "$rev_xmin $rev_xmax  $rev_parhets  $rev_xmin_norm $rev_xmax_norm ",
      "\n";
  }
}
