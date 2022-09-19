#!/usr/bin/perl -w
use strict;

# read in a cluster file (output of agmr_cluster ) and a
# pedigree file (tab separated, accession, Fparent, Mparent ids in last 3 cols of each line)
# output to stdout a pedigree file with the parental ids replaced with the representative of the cluster

my $cluster_file = shift;
my $pedigree_file = shift;

# store cluster info in hash
my %clusterids_repid = (); # key: an id in a cluster, value: the representative id for that cluster
my %repids = (); # keys are ids of the representative ids of the clusters (with size >= 2)
open my $fhin, "<", "$cluster_file";
while(my $line = <$fhin>){
    next if($line =~ /^\s*#/);
    my @cols = split(" ", $line);
    my $the_repid = $cols[4];
    $repids{$the_repid} = 1;
    for my $id (@cols[4..$#cols]){
     # print STDERR "$id $the_repid\n";
	$clusterids_repid{$id} = $the_repid;
    }
}
close $fhin;
print STDERR "# n clusters: ", scalar keys %repids, "  n ids: ", scalar keys %clusterids_repid, "\n";

open $fhin, "<", "$pedigree_file";
# my $line1 = <$fhin>;
# print "xxx: ", $line1;
while(my $line = <$fhin>){
  next if($line =~ /^\s*#/);
  $line =~ s/\s+$//; # remove newline, any other whitespace at end.
    my @cols = split(" ", $line);
    # replace the last 2 elements of @cols with cluster representatives, if belong to cluster
    #   print STDERR $cols[-1], "  ", $cols[-2], "  ";
    # my $mparent = $cols[-1];
    # my $fparent = $cols[-2];
    # if($mparent =~ /\s/  or  $fparent =~ /\s/){
    #   print STDERR "[$fparent] [$mparent] \n";
    # }
    $cols[-1] = $clusterids_repid{$cols[-1]} // $cols[-1];
    $cols[-2] = $clusterids_repid{$cols[-2]} // $cols[-2];
 #   print STDERR $cols[-1], "\t", $cols[-2], "\n";
    print join(" ", @cols), "\n";
}
close $fhin;
