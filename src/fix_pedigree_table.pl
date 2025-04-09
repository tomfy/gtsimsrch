#!/usr/bin/perl -w
use strict;
use List::Util qw(min max sum);

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
    next if(/^\s*#|^\s*$/);
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

my $n_complete_pedigrees_1 = 0;
my $n_complete_pedigrees_2 = 0;
my %u_progeny_ids = ();
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
  next if($A eq 'NA'  or  $F eq 'NA'  or $M eq 'NA');
  if(!exists $u_progeny_ids{$A}){
    $u_progeny_ids{$A} = [];
  }
  $n_complete_pedigrees_1++;
  push @{$u_progeny_ids{$A}}, "$F  $M";
  print "$A  $F  $M\n";
}


while(my($progid, $peds) = each %u_progeny_ids){
  $n_complete_pedigrees_2 += scalar @$peds;
  print STDERR "$progid\n  ", join("\n  ", @$peds), "\n" if(scalar @$peds > 1);
}
print STDERR "n distinct uniq. accessions with pedigrees: ", scalar keys %u_progeny_ids, "\n";
print STDERR "n pedigrees:  $n_complete_pedigrees_1  $n_complete_pedigrees_2 \n";
