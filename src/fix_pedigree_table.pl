#!/usr/bin/perl -w
use strict;

# take the last 3 fields, and print them as first three
# AND, if a cluster file is given
# replace the accessions ids with those of the
# representative accessions of the duplicate cluster

# usage:  fix_pedigree_table.pl  <cluster_filename>  <  <pedigree_table_filename>

my $cluster_file_name = shift // undef;
my %id_repid = ();
if(defined $cluster_file_name){
  open my $fhcluster, "<", "$cluster_file_name";
  while(<$fhcluster>){
    next if(/^\s*#/);
    my @cluster_ids = split(" ", $_); # 9 cols of other stuff, then ids of cluster members
    my $cluster_size = $cluster_ids[0];
    @cluster_ids = @cluster_ids[9..$#cluster_ids];
    die if(scalar @cluster_ids != $cluster_size);
    my $rep_id = $cluster_ids[0];
    for my $an_id (@cluster_ids){
      $id_repid{$an_id} = $rep_id;
    }
  }
}

my $x = <>; # skip first line
while(<>){
  next if(/^\s*#/); # skip comments
  s/\s+$//;
  my @fields = split(" ", $_);
  my ($A, $F, $M) = @fields[-3,-2,-1];
  # if ids are members of clusters, replace with representative accessions id.
  $A = $id_repid{$A} // $A;
  $F = $id_repid{$F} // $F;
  $M = $id_repid{$M} // $M;
  
  print "$A  $F  $M\n";
}
