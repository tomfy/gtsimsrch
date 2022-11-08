#!/usr/bin/perl -w
use strict;

# runs simsearch, then agmr_cluster, and outputs a file
# with the same format as simsearch input, but with just one
# line representing each cluster.

my $input_dosages_filename = shift;
my $max_acc_missing_data_fraction = shift // 0.5;
open my $fhin, "<", "$input_dosages_filename" or die "couldn't open $input_dosages_filename for reading.\n";
my $cleaned_dosages_filename = $input_dosages_filename . "_cleaned";
open my $fhout, ">", "$cleaned_dosages_filename";
my $n_markers = undef;
my $n_bad_accessions = 0;
print STDERR "# Removing accession with too much missing data.\n";
while(my $line_in = <$fhin>){
  if($line_in =~ /^\s*#/){
   # print STDERR "print comment line.\n";
    print $fhout $line_in;
  }elsif($line_in =~ /^MARKER/){
  #  print STDERR "print marker ids line.\n";
    print $fhout $line_in;
    my @markers = split(" ", $line_in);
    $n_markers = scalar @markers  - 1;
  }else{
    die "File lacks line with MARKER and marker ids\n" if(! defined $n_markers);
    my $n_bad = () = $line_in =~ /\sNA/g;
    if($n_bad/$n_markers <= $max_acc_missing_data_fraction){
      print $fhout $line_in;
    }else{
      $n_bad_accessions++;
    }
  }
}
close $fhout;
print STDERR "# $n_bad_accessions accessions eliminated due to excessive missing data.\n";
print "# $n_bad_accessions accessions eliminated due to excessive missing data (>" ,
  int($max_acc_missing_data_fraction*100 + 0.5), "\%)\n";

my $simsearch_command = "simsearch  -i $cleaned_dosages_filename";

system "$simsearch_command";

# run agmr_cluster to get clusters of near-identical genotypes
system "agmr_cluster < simsearch.out > agmr_cluster.out";

# store individual lines of dosage file in hash
open my $fh_dosages, "<", "$cleaned_dosages_filename";

my %id_line = ();
my %id_mdcount = ();
my $first_line = <$fh_dosages>;
while(my $line = <$fh_dosages>){
    next if($line =~ /^\s*#/); # skip comments
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

while (my $line = <$fh_clusters>) {
  next if($line =~ /^\s*#/);
  my @cols = split(" ", $line);
  my $min_md_id = $cols[4];
  my $min_md_count = $id_mdcount{$min_md_id};
 # if (1) { # uniquified file has accession with the least missing data.
    for my $a_cluster_id (@cols[4..$#cols]) {
      $clusterids{$a_cluster_id} = 1; # all accessions in clusters get stored in @clusterids
      if (0) {
	if ($id_mdcount{$a_cluster_id} < $min_md_count) { # find the member of cluster with the least missing data
	  $min_md_count =  $id_mdcount{$a_cluster_id};
	  $min_md_id = $a_cluster_id;
	}
      }
    }
    print $id_line{$min_md_id};
  # } else { # uniquified file just has first-occurring accession in cluster 
  #   print $id_line{$cols[4]};
  # }
}

# print the lines of the singletons (accessions in clusters of size 1)
my @sorted_ids = sort {$a cmp $b} keys %id_line; # sort keys and output in this order for repeatability
#while (my($an_id, $a_line) = each %id_line) {
for my $an_id (@sorted_ids){
  if (! exists $clusterids{$an_id}) {
    print $id_line{$an_id};
  }
}

