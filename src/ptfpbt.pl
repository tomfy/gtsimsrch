#!/usr/bin/perl -w
use strict;

my $ptfile = shift;
my $btfile = shift;
my $maxD = shift // 0.03;

open my $fhpt, "<", "$ptfile";
open my $fhbt, "<", "$btfile";

my %accid_ptline = ();
my %accid_ptppd = ();
while(my $line = <$fhpt>){
  my @cols = split(" ", $line);
  my $acc = $cols[0];
  my $D = $cols[17];
  if($D <= $maxD){
    $accid_ptline{$acc} = $line;
    my ($p1, $p2) = @cols[2,3];
    my $ppd = order_pair($p1, $p2) . " $D";
    $accid_ptppd{$acc} = $ppd;
  }

}

while(my $line = <$fhbt>){
  my @cols = split(" ", $line);
  my $accid = $cols[0];
  my $ptppd = $accid_ptppd{$accid} // undef;
  if(defined $ptppd){
  for(my ($offset, $i) = (2, 0); $offset+17 < scalar @cols; $offset += 23, $i++){
    my $ppd = get_p1p2D(\@cols, $offset);
    my $best_ppd = $ppd if($offset == 2);
  }
}else{

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
