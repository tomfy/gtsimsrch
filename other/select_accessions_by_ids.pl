#!/usr/bin/perl/ -w
use strict;

#

my $idfile = shift;
my $gtfile = shift;
my $idcol = shift // 1;
$idcol--;

my %id = ();
open my $fhin, "<", $idfile;
while(<$fhin>){
  next if(/^\s*#/);
  my @cols = split(" ", $_);
  if(scalar @cols >= $idcol){
    $id{$cols[$idcol]}++;
  }
}
close $fhin;

# while(my ($id, $c) = each %id){
#   print STDERR "[$id] $c \n";
# }
# getc();

open $fhin, "<", $gtfile;
while(<$fhin>){ # skip initial comment lines
  next if(/^\s*#/);
  print if(/^\s*MARKER/); # print 1st non-comment line (with marker ids)
  if(/^\s*CHROMOSOME/){
    print;
    last;
  }
}

while(<$fhin>){
  next if(/^\s*#/);
  if(/^\s*(\S+)/){
   #print "{$1}\n";
    if(exists $id{$1}){
      print;
      #print STDERR "$1\n";
      #getc();
    }
  }
}
