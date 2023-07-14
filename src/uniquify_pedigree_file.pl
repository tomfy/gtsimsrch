#!/usr/bin/perl -w
use strict;
use Getopt::Long;

# read in a cluster file (output of agmr_cluster ) and a
# pedigree file (whitespace separated, accession, Fparent, Mparent ids in last 3 cols of each line)
# output to stdout a pedigree file with the parental ids replaced with the representative of the cluster

my $cluster_file = undef;
my $pedigree_file = undef;

GetOptions(
	   'cluster_file=s' => \$cluster_file,
	   'pedigree_file=s' => \$pedigree_file,
	  );

# store cluster info in hash
my %clusterid_repid = (); # key: an id in a cluster, value: the representative id for that cluster
my %repids = (); # keys are ids of the representative ids of the clusters (with size >= 2)
open my $fhin, "<", "$cluster_file";
while (my $line = <$fhin>) {
  next if($line =~ /^\s*#/);
  my @cols = split(" ", $line);
  my $the_repid = $cols[9];
  $repids{$the_repid} = 1;
  for my $id (@cols[9..$#cols]) {
    # print STDERR "$id $the_repid\n";
    $clusterid_repid{$id} = $the_repid;
  }
}
close $fhin;
print STDERR "# n clusters: ", scalar keys %repids, "  n ids: ", scalar keys %clusterid_repid, "\n";

open $fhin, "<", "$pedigree_file";
my $x = <$fhin>;		# first line - throw away.
# my $line1 = <$fhin>;
# print "xxx: ", $line1;
while (my $line = <$fhin>) {
  next if($line =~ /^\s*#/);
  $line =~ s/\s+$//;	# remove newline, any other whitespace at end.
  my @cols = split(" ", $line);
  my $progeny_id = $cols[-3];

  # if ($progeny_id ne 'NA') {
  #   if ((exists $clusterid_repid{$progeny_id}) and ($clusterid_repid{$progeny_id} ne $progeny_id)) { # if progeny_id belongs to a cluster, only output if it is the representative accession.
  #     print STDERR "$progeny_id is non-rep member of cluster with rep ", $clusterid_repid{$progeny_id}, "\n";
  #   } else {
  #     my $Fpar_id = $clusterid_repid{$cols[-2]} // $cols[-2];
  #     my $Mpar_id = $clusterid_repid{$cols[-1]} // $cols[-1];
  #     print "$progeny_id $Fpar_id $Mpar_id\n";
  #   }
  # }
  if($progeny_id ne 'NA'){ # output, replacing ids by ids of the cluster representatives if in a cluster 
    $progeny_id = $clusterid_repid{$progeny_id} // $progeny_id;
     my $Fpar_id = $clusterid_repid{$cols[-2]} // $cols[-2];
    my $Mpar_id = $clusterid_repid{$cols[-1]} // $cols[-1];
    print "$progeny_id  $Fpar_id  $Mpar_id\n";
  }
}
close $fhin;
