#!/usr/bin/perl -w
use strict;

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
	$clusterids_repid{$id} = $the_repid;
    }
}
close $fhin;
print STDERR "# n clusters: ", scalar keys %repids, "  n ids: ", scalar keys %clusterids_repid, "\n";

open $fhin, "<", "$pedigree_file";
my $line1 = <$fhin>;
print $line1;
while(my $line = <$fhin>){
    chomp $line;
    my @cols = split("\t", $line);
    # replace the last 2 elements of @cols with cluster representatives, if belong to cluster
 #   print STDERR $cols[-1], "  ", $cols[-2], "  ";
    $cols[-1] = $clusterids_repid{$cols[-1]} // $cols[-1];
    $cols[-2] = $clusterids_repid{$cols[-2]} // $cols[-2];
 #   print STDERR $cols[-1], "\t", $cols[-2], "\n";
    print join("\t", @cols), "\n";
}
close $fhin;
