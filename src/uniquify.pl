#!/usr/bin/perl -w
use strict;

# run simsearch, then cluster output, and output a file
# with the same format as simsearch input, but with just one
# line representing each cluster.

my $input_dosages_file = shift;

my $simsearch_command = "simsearch  -i $input_dosages_file";

system "$simsearch_command";

system "agmr_cluster < simsearch.out > agmr_cluster.out";

open my $fh_dosages, "<", "$input_dosages_file";

my %id_line = ();
my $first_line = <$fh_dosages>;
while(my $line = <$fh_dosages>){
    if($line =~ /^\s*(\S+)/){
	$id_line{$1} = $line;
    }
}
close($fh_dosages);

open my $fh_clusters, "<", "agmr_cluster.out";

my %clusterids = ();
print $first_line;
while(my $line = <$fh_clusters>){
    next if($line =~ /^\s*#/);
    my @cols = split(" ", $line);
    my $representative_id = $cols[4];
    print $id_line{$representative_id} // "####\n";
    for my $a_cluster_id (@cols[4..$#cols]){
	$clusterids{$a_cluster_id} = 1; # value doesn't matter, key just needs to exist in hash.
    }
}

while(my($an_id, $a_line) = each %id_line){
    if(! exists $clusterids{$an_id}){
	print $a_line;
    }
}

