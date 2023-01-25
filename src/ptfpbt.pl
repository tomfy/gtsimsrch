#!/usr/bin/perl -w
use strict;
use Scalar::Util qw(looks_like_number);

my $ptfile = shift;
my $btfile = shift;
my $maxD = shift // 0.03;

open my $fhpt, "<", "$ptfile";
open my $fhbt, "<", "$btfile";

my %accid_ptline = ();
my %accid_ptppd = ();
# store pedigree_test information
while (my $line = <$fhpt>) {
  my @cols = split(" ", $line);
  my $acc = $cols[0];
  my $D = $cols[17];
  if (looks_like_number($D)  and  $D <= $maxD) {
    $accid_ptline{$acc} = $line;
    my ($p1, $p2) = @cols[2,3];
    my $ppd = order_pair($p1, $p2) . " $D"; # Fparent id, Mparent id, D
    $accid_ptppd{$acc} = $ppd;
  }
}

# store find_parents/best_triples information
while (my $line = <$fhbt>) {
  my @cols = split(" ", $line);
  my $accid = $cols[0];
  my $ptppd = $accid_ptppd{$accid} // undef;
  #my ($ptp1, $ptp2, $ptD) = split(" ", $ptppd);
  my $best_ppd = "XXX";
  my $next_best_ppd = "xxx";
  my $pedigree_rank = -1;
  if (defined $ptppd) { # if pedigree_test had result for this accession ...
    for (my ($offset, $i) = (2, 0); $offset+17 < scalar @cols; $offset += 23, $i++) {
      my $ppd = get_p1p2Dstd(\@cols, $offset);
     # print STDERR "***   $accid  $ppd \n";
      if($i == 0){
	$best_ppd = $ppd;
      }elsif($i == 1){
	$next_best_ppd = $ppd;
      }
      my $ptpp = ($ptppd =~ /^(\S+\s+\S+)/)? $1 : 'x x';
      my $fppp = ($ppd =~ /^(\S+\s+\S+)/)? $1 : 'x x';
      if($ptpp eq $fppp){ # just check (ordered) parental ids are the same
	$pedigree_rank = $i;
	last if($i >= 1);
      }
    }
    print "$accid    $ptppd   $best_ppd   $next_best_ppd   $pedigree_rank \n"; 
  } else {
    print STDERR "No pedigree for accession: $accid \n";
  }
}


sub get_p1p2Dstd{
  my $cols = shift;
  my $offset = shift;

  my $p1 = $cols->[$offset+1];
  my $p2 = $cols->[$offset+2];

  my $parent_pair = order_pair($p1, $p2);
  my $D = $cols->[$offset+17];
  return "$parent_pair $D";
}

sub order_pair{
  my $p1 = shift;
  my $p2 = shift;
  return ($p1 le $p2)? "$p1 $p2" : "$p2 $p1";
}
