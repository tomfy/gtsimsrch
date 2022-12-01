#!/usr/bin/perl -w
use strict;

# runs simsearch, then agmr_cluster, and outputs a file
# with the same format as simsearch input, but with just one
# line representing each cluster.


my $input_dosages_filename = shift;
my $do_remove_bad_accessions = shift // 1;
my $max_acc_missing_data_fraction = shift // 0.5;
my $simsearch_command;
open my $fhin, "<", "$input_dosages_filename" or die "couldn't open $input_dosages_filename for reading.\n";
my $cleaned_dosages_filename = $input_dosages_filename . "_cleaned";

# remove accessions with excessive missing data
if ($do_remove_bad_accessions) {
  open my $fhout, ">", "$cleaned_dosages_filename";
  my $n_markers = undef;
  my $n_bad_accessions = 0;
  print STDERR "# Removing accessions with too much missing data.\n";
  while (my $line_in = <$fhin>) {
    if ($line_in =~ /^\s*#/) {
      # print STDERR "print comment line.\n";
      print $fhout $line_in;
    } elsif ($line_in =~ /^MARKER/) {
      #  print STDERR "print marker ids line.\n";
      print $fhout $line_in;
      my @markers = split(" ", $line_in);
      $n_markers = scalar @markers  - 1;
    } else {
      die "File lacks line with MARKER and marker ids\n" if(! defined $n_markers);
      my $n_bad = () = $line_in =~ /\sNA/g;
      if ($n_bad/$n_markers <= $max_acc_missing_data_fraction) {
	print $fhout $line_in;
      } else {
	$n_bad_accessions++;
      }
    }
  }
  close $fhout;
  print STDERR "# $n_bad_accessions accessions eliminated due to excessive missing data.\n";
  print "# $n_bad_accessions accessions eliminated due to excessive missing data (>" ,
    int($max_acc_missing_data_fraction*100 + 0.5), "\%)\n";

  $simsearch_command = "simsearch  -i $cleaned_dosages_filename";

} else {
  $simsearch_command = "simsearch -i $input_dosages_filename";
}

system "$simsearch_command";

# run agmr_cluster to get clusters of near-identical genotypes
system "agmr_cluster < simsearch.out > agmr_cluster.out";

open my $fh_clusters, "<", "agmr_cluster.out";
my %clusterids = ();

while (my $line = <$fh_clusters>) { # each line is one cluster
  next if($line =~ /^\s*#/);
  my @cols = split(" ", $line);
  my $cluster_size = shift @cols;
  my $min_d = shift @cols;
  my $max_d = shift @cols;
  my $n_bad = shift @cols;
  my $rep_id = $cols[0];     # id of the representative of the cluster
  for my $a_cluster_id (@cols) {
    $clusterids{$a_cluster_id} = 1; # all accessions in clusters get stored in @clusterids
  }
  print STDERR "done storing cluster of size $cluster_size \n";
  #  my $elected_gts = vote(\@cols, \%id_gts);
  #  print STDERR "done with cluster vote \n";
  #  print "$rep_id  ", join(" ", @$elected_gts), "\n";

}
close $fh_clusters;


print STDERR "before storing dosages\n";

# store individual lines of dosage file in hash
open my $fh_dosages, "<", "$cleaned_dosages_filename";

# store ids and genotypes of clusters (size >= 2), and output singletons.
my %id_gts  = ();  # key: ids; value: array ref of dosages (0,1,2,NA) 
my $first_line = <$fh_dosages>;
print $first_line;
while (my $line = <$fh_dosages>) {
  next if($line =~ /^\s*#/);	# skip comments
  my @cols = split(" ", $line);
  my $id = shift @cols;
  if (exists $clusterids{$id}) {
    $id_gts{$id} = \@cols;
  } else {
    print $line;
  }
}
close($fh_dosages);

print STDERR "after storing dosages\n";
# $x = getc();

open  $fh_clusters, "<", "agmr_cluster.out";
while (my $line = <$fh_clusters>) { # each line is one cluster
  next if($line =~ /^\s*#/);
  my @cols = split(" ", $line);
  my $cluster_size = shift @cols;
  my $min_d = shift @cols;
  my $max_d = shift @cols;
  my $n_bad = shift @cols;
  my $rep_id = $cols[0];     # id of the representative of the cluster

  #  print STDERR "done storing cluster of size $cluster_size \n";
  print STDERR "before vote\n";
  my $elected_gts = vote(\@cols, \%id_gts);
  print STDERR "size of elected_gts: ", scalar @$elected_gts, "\n";
  #  print STDERR "done with cluster vote \n";
  print "$rep_id  ", join(" ", @$elected_gts), "\n";

}

# while (my $line = <$fh_clusters>) { # each line is one cluster
#   next if($line =~ /^\s*#/);
#   my @cols = split(" ", $line);
#   my $min_md_id = $cols[4];
#   my $min_md_count = $id_mdcount{$min_md_id};
#  # if (1) { # uniquified file has accession with the least missing data.
#     for my $a_cluster_id (@cols[4..$#cols]) {
#       $clusterids{$a_cluster_id} = 1; # all accessions in clusters get stored in @clusterids
#       if (0) {
# 	if ($id_mdcount{$a_cluster_id} < $min_md_count) { # find the member of cluster with the least missing data
# 	  $min_md_count =  $id_mdcount{$a_cluster_id};
# 	  $min_md_id = $a_cluster_id;
# 	}
#       }
#     }
#     print $id_line{$min_md_id};
#   # } else { # uniquified file just has first-occurring accession in cluster 
#   #   print $id_line{$cols[4]};
#   # }
# }

# print the lines of the singletons (accessions in clusters of size 1)
# my @sorted_ids = sort {$a cmp $b} keys %id_gts; # sort keys and output in this order for repeatability
# #while (my($an_id, $a_line) = each %id_line) {
# for my $an_id (@sorted_ids) {
#   if (! exists $clusterids{$an_id}) {
#     print "$an_id  ", join(" ", $id_gts{$an_id}), "\n";
#   }
# }

# ************************************************

sub vote{
  my $cluster_ids = shift;	# array ref holding the ids of the accessions in the cluster.
  my $id_dosages = shift;	# hash ref. key: id, value: array ref of dosages for this id.
print STDERR "in vote. new cluster:  \n";
  my $first_id = $cluster_ids->[0];
  my $first_dosages = $id_dosages->{$first_id};
  my $n_markers = scalar @$first_dosages;
  my $cluster_size = scalar @$cluster_ids;
  my @marker_votes = ();;
  for(1..$n_markers){
push @marker_votes, [0, 0, 0, 0];
  }
  for my $an_id (@$cluster_ids) {
    my $dosages = $id_dosages->{$an_id}; # array ref of dosages for $an_id
    while (my($i, $d) = each @$dosages) {
      $d = 3 if($d eq 'NA');
      $marker_votes[$i]->[$d]++;
    }
 #   print STDERR "  $an_id   ";
    # for my $jjj (0..9){
    #   my $vs = $marker_votes[$jjj];
    #   print STDERR join(" ", @$vs), ";  ";
    # }print STDERR "\n";
  }
#  my $x = getc();
  print STDERR "in vote. after storing votes of all markers, all accessions in cluster. \n";
  my @elected_dosages = (0) x $n_markers;
  while (my($i, $mv) = each @marker_votes) {
    my $e = 'NA';
 #   print STDERR "    $i  ", join(" ", @$mv), "\n";
    for my $j (0..2) {
      my $vote = $mv->[$j]; # the number of votes for dosage = $j for marker $i
      if (2*$vote > $cluster_size) { # must have > 50% for one dosage, or 'NA'
	$e = $j;
	last;
      }
    }
    $elected_dosages[$i] = $e;
  }
  return \@elected_dosages;
}
