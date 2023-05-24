#!/usr/bin/perl -w
use strict;

my $file1 = shift;
my $file2 = shift;

my %idpair_agmr;
open my $fh1, "<", "$file1";
while(my $line = <$fh1>){
  next if($line =~ /^\s*#/);
  my @cols = split(" ", $line);
  my ($id1, $id2, $agmr) = @cols[0, 1, 5];
 # print "$id1  $id2  $agmr\n";
  my $idpair = ($id1 le $id2)? "$id1 $id2" : "$id2 $id1";
  $idpair_agmr{$idpair} = $agmr;
}


open my $fh2, "<", "$file2";
while(my $line = <$fh2>){
  next if($line =~ /^\s*#/);
  my @cols = split(" ", $line);
  my ($id1, $id2, $agmr) = @cols[0, 1, 5];
  my $idpair = ($id1 le $id2)? "$id1 $id2" : "$id2 $id1";
  if(exists $idpair_agmr{$idpair}){
    print "$idpair  ", $idpair_agmr{$idpair}, "  $agmr \n";
  }
}
