#!/usr/bin/perl -w
use strict;

# read in pedigree file (3 columns: accession Fparent Mparent)
# find accession triples which are (according to pedigree file)  accession child parent, i.e. col3 is grandparent of col2, col1 is parent of col2 and child of col 3
my %A_F = ();
my %A_M = ();
while(<>){
  next if(/^\s*#/);
  my @cols = split(" ", $_);
  my ($A, $F, $M) = @cols[0,1,2];
  if($F ne 'NA'){
    $A_F{$A} = $F;
  }
  if($M ne 'NA'){
    $A_M{$A} = $M;
  }
}

print "#  accession  parent  grandparent \n";
while(my($a, $f) = each %A_F){
  my $gp = $A_F{$f} // undef;
  if(defined $gp){
    print "$a $f  $gp   maternal grandma\n";
  }
  $gp = $A_M{$f} // undef;
  if(defined $gp){
    print "$a  $f  $gp   maternal grandpa\n";
  }
}

while(my($a, $m) = each %A_M){
  my $gp = $A_F{$m} // undef;
  if(defined $gp){
    print "$a  $m  $gp   paternal grandma\n";
  }
  $gp = $A_M{$m} // undef;
  if(defined $gp){
    print "$a  $m  $gp   paternal grandpa\n";
  }
}
  
