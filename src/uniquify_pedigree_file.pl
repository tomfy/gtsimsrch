#!/usr/bin/perl -w
use strict;

# read in a cluster file (output of agmr_cluster ) and a
# pedigree file (whitespace separated, accession, Fparent, Mparent ids in last 3 cols of each line)
# output to stdout a pedigree file with the parental ids replaced with the representative of the cluster

my $cluster_file = shift;
my $pedigree_file = shift;

# store cluster info in hash
my %clusterids_repid = (); # key: an id in a cluster, value: the representative id for that cluster
my %repids = (); # keys are ids of the representative ids of the clusters (with size >= 2)
open my $fhin, "<", "$cluster_file";
while (my $line = <$fhin>) {
  next if($line =~ /^\s*#/);
  my @cols = split(" ", $line);
  my $the_repid = $cols[4];
  $repids{$the_repid} = 1;
  for my $id (@cols[4..$#cols]) {
    # print STDERR "$id $the_repid\n";
    $clusterids_repid{$id} = $the_repid;
  }
}
close $fhin;
print STDERR "# n clusters: ", scalar keys %repids, "  n ids: ", scalar keys %clusterids_repid, "\n";

open $fhin, "<", "$pedigree_file";
my $x = <$fhin>;		# first line - throw away.
# my $line1 = <$fhin>;
# print "xxx: ", $line1;
while (my $line = <$fhin>) {
  next if($line =~ /^\s*#/);
  $line =~ s/\s+$//;	# remove newline, any other whitespace at end.
  my @cols = split(" ", $line);
  my $prog_id = $cols[-3];

  # if ($prog_id ne 'NA') {
  #   if ((exists $clusterids_repid{$prog_id}) and ($clusterids_repid{$prog_id} ne $prog_id)) { # if prog_id belongs to a cluster, only output if it is the representative accession.
  #     print STDERR "$prog_id is non-rep member of cluster with rep ", $clusterids_repid{$prog_id}, "\n";
  #   } else {
  #     my $Fpar_id = $clusterids_repid{$cols[-2]} // $cols[-2];
  #     my $Mpar_id = $clusterids_repid{$cols[-1]} // $cols[-1];
  #     print "$prog_id $Fpar_id $Mpar_id\n";
  #   }
  # }
  if($prog_id ne 'NA'){
    $prog_id = $clusterids_repid{$prog_id} // $prog_id;
     my $Fpar_id = $clusterids_repid{$cols[-2]} // $cols[-2];
    my $Mpar_id = $clusterids_repid{$cols[-1]} // $cols[-1];
    print "$prog_id  $Fpar_id  $Mpar_id\n";
  }
}
close $fhin;
