#!/usr/bin/perl -w
use strict;

my $file1 = shift;
my $col1 = shift;
my $file2 = shift;
my $col2 = shift;
$col1--;
$col2--;

my %idpair_agmr1;
open my $fh1, "<", "$file1";
while(my $line = <$fh1>){
  next if($line =~ /^\s*#/);
  my @cols = split(" ", $line);
  my ($id1, $id2, $agmr) = @cols[0, 1, $col1];
 # print "$id1  $id2  $agmr\n";
  my $idpair = ($id1 le $id2)? "$id1 $id2" : "$id2 $id1";
  $idpair_agmr1{$idpair} = $agmr;
}

my $f2notf1 = ""; my $count2not1 = 0;
my $f1notf2 = ""; my $count1not2 = 0;
open my $fh2, "<", "$file2";
my %idpair_agmr2 = ();
while(my $line = <$fh2>){
  next if($line =~ /^\s*#/);
  my @cols = split(" ", $line);
  my ($id1, $id2, $agmr) = @cols[0, 1, $col2];
  my $idpair = ($id1 le $id2)? "$id1 $id2" : "$id2 $id1";
  $idpair_agmr2{$idpair} = $agmr;
  if(exists $idpair_agmr1{$idpair}){
    print "$idpair  ", $idpair_agmr1{$idpair}, "  $agmr \n";
  }else{ # in file2, but not file1
    $f2notf1 .= "$idpair $agmr \n";
    $count2not1++;
  }
}

while(my($idp, $a) = each %idpair_agmr1){
  if(!exists $idpair_agmr2{$idp}){
    $f1notf2 .= "$idp $a \n";
    $count1not2++;
  }
}

print "## $count1not2  $count2not1\n";

print STDERR $f2notf1;

