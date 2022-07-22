#!/usr/bin/perl -w
use strict;

# run simsearch, then cluster output, and output a file
# with the same format as simsearch input, but with just one
# line representing each cluster.

my $input_dosages_file = shift;

my $simsearch_command = "simsearch  -i $input_dosages_file";

system "$simsearch_command";

# run agmr_cluster to get clusters of near-identical genotypes
system "agmr_cluster < simsearch.out > agmr_cluster.out";

# store individual lines of dosage file in hash
open my $fh_dosages, "<", "$input_dosages_file";

my %id_line = ();
my %id_mdcount = ();
my $first_line = <$fh_dosages>;
while(my $line = <$fh_dosages>){
    next if($line =~ /^\s*#/);
    if($line =~ /^\s*(\S+)/){
	my $id = $1;
	$id_line{$id} = $line;
	$line =~ s/^\s*\S+//;
	my $md_count = () = ($line =~ /NA/g);
	$id_mdcount{$id} = $md_count;
#	print $md_count, "  ", $line;
    }
}
close($fh_dosages);


open my $fh_clusters, "<", "agmr_cluster.out";
my %clusterids = ();
print $first_line;

while(my $line = <$fh_clusters>){
    next if($line =~ /^\s*#/);
    my @cols = split(" ", $line);
    my $min_md_id = $cols[4];
    my $min_md_count = $id_mdcount{$min_md_id};
    if(0){ # uniquified file has accession with the least missing data.
    for my $a_cluster_id (@cols[4..$#cols]){
	$clusterids{$a_cluster_id} = 1; # value doesn't matter, key just needs to exist in hash.
	if($id_mdcount{$a_cluster_id} < $min_md_count){ # find the member of cluster with the least missing data
	    $min_md_count =  $id_mdcount{$a_cluster_id};
	    $min_md_id = $a_cluster_id;
	}
    }
    print $id_line{$min_md_id};
    }lse{ # uniquified file just has first-occurring accession in cluster 
	print $id_line{$cols[4]};
    }
}

# print the lines of the singletons (accessions in clusters of size 1)
while(my($an_id, $a_line) = each %id_line){
    if(! exists $clusterids{$an_id}){
	print $a_line;
    }
}

