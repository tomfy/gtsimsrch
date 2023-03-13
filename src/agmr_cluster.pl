#!/usr/bin/perl -w
use strict;
use Graph::Undirected;
use Getopt::Long;
use File::Basename 'dirname';

use Cwd 'abs_path';
my ( $bindir, $libdir );
BEGIN {     # this has to go in Begin block so happens at compile time
  $bindir =
    dirname( abs_path(__FILE__) ) ; # the directory containing this script
  $libdir = $bindir . '/../moose/lib';
  $libdir = abs_path($libdir);	# collapses the bin/../lib to just lib
}
use lib $libdir;
use Cluster1d;

# read suplicatesearch output and find clusters of accessions
# with near-identical genotypes
# specifically, create a graph, and for any pairs in duplicatesearch output with
# agmr <= $cluster_max_agmr, make an edge in the graph between the accessions of the pair.
# the clusters are then the connected components of the graph.
# i.e. two accession belong to the same cluster iff you can get from
# one to the other by traversing edges of the graph.
# one line of output for each cluster, with fields
# cluster size, min agmr, max agmr, 'nbad (should be 0)', the ids of accessions in cluster.
# e.g.:    3  0.0104 0.0142 0  TMEB778:250254008  TMEB778:250304613  TMEB778:250597946
# the duplicatesearch output should have all the agmrs for all the pairs in each cluster,
# and in that case the 4th field ('nbad') will be 0. If nbad > 0,
# consider rerunning this script with a smaller $max_agmr,
# or rerun duplicatesearch with a larger value of 'max_estimated_agmr' (-e option)

# usage example:
# agmr_cluster.pl -in duplicatesearch.out  -out  agmr_clusters  [-cluster_max 0.08]

# runs duplicatesearch, then agmr_cluster, and outputs a file
# with the same format as duplicatesearch input, but now with just one
# line representing each cluster.

my $input_agmr_filename = undef;
my $cluster_max_agmr = 'auto'; # construct graph with edges between pairs of accessions iff their agmr is <= this.
my $output_cluster_filename = "agmr_cluster.out";
my $pow = 1;			# 'log';
my $minx = 0.001;
my $column = 6; # the agmrs to use for cluster are found in this column (unit-based)

GetOptions(
	   'input_file=s' => \$input_agmr_filename, # file with id1 id2 x xx agmr_est agmr (duplicatesearch output)
	   'output_file=s' => \$output_cluster_filename,
	   'cluster_max_agmr=s' => \$cluster_max_agmr, # cluster using graph with edges for pairs with agmr < this.
	   'pow=s' => \$pow,
	   'column=i' => \$column,
	  );
print STDERR "Clustering based on column: $column\n";
$column--;
if (!defined $input_agmr_filename) {
  print STDERR "Basic usage example: \n", "agmr_cluster -in duplicatesearch.out  -out acluster.out \n";
  print STDERR "by default agmr_cluster will attempt to automatically decide the max agmr between duplicates.\n",
    " but you can specify it with the cluster option, e.g.  -cluster 0.06 \n";
  exit;
}

# store all agmrs in the input file in hash %edge_weight
my %edge_weight = (); # keys are ordered pairs of accession ids representing graph edges, values are agmrs
my %id_closeidas = (); # keys are ids, values ids and agmrs of other close accessions. id1:["$id2 $agmr12", "$id3 $agmr13", ...]
#my %id_agmrs = (); # keys are ids, values array ref of agmrs found between that id and others
open my $fhin, "<", "$input_agmr_filename" or die "Couldn't open $input_agmr_filename for reading.\n";
while (my $line = <$fhin>) {
  next if($line =~ /^\s*#/);
  my @cols = split(" ", $line);
  my ($id1, $id2, $usable_chunks, $match_chunks, $est_agmr) = @cols[0..4];
  my $agmr = $cols[$column];
  my $edge_verts = ($id1 lt $id2)? "$id1 $id2" : "$id2 $id1"; # order the pair of ids
  if (!exists $id_closeidas{$id1}) {
    $id_closeidas{$id1} = ["$id2 $agmr"];
  } else {
    push @{$id_closeidas{$id1}}, "$id2 $agmr";
  }
  # if(!exists $id_agmrs{$id1}){
  #   $id_agmrs{$id1} = [$agmr];
  # }else{
  #   push @{$id_agmrs{$id1}}, $agmr;
  # }
  $edge_weight{$edge_verts} = $agmr;
}
close $fhin;

my @agmrs = values %edge_weight;
if ($cluster_max_agmr eq 'auto') {
 
  my $cluster1d_obj = Cluster1d->new({label => '', xs => \@agmrs, pow => $pow, minx => $minx});
  print  "before two_cluster \n";
  print  "#  pow: $pow  min: $minx \n";
  my ($Hopt, $maxQ, $Hmid_half_max) = $cluster1d_obj->two_cluster();
  # my ($Hoptx, $maxQx) = (0, 0); # $cluster1d_obj->two_cluster_x();
  print "# after two_cluster \n";
  #  print "# before k-means/kde clustering \n";
  #  my ($n_pts, $km_n_L, $km_n_R, $km_h_opt, $q, $kde_n_L, $kde_n_R, $kde_h_opt, $kde_q) = $cluster1d_obj->one_d_2cluster();
  #  printf( STDERR "# clustering %5d points;  k-means: %5d below  %8.6f and  %5d above; q: %6.4f.  kde: %5d below  %8.6f  and %5d above; kde_q: %6.4f   Hopt: %6.4f  maxQ: %6.4f.  Hmhmx: %6.4f \n",
  #       $n_pts, $km_n_L, $km_h_opt, $km_n_R, $q, $kde_n_L, $kde_h_opt, $kde_n_R, $kde_q, $Hopt, $maxQ, $Hmid_half_max);
  printf(STDERR  "Hopt: %6.4f  maxQ: %6.4f.  Hmhmx: %6.4f \n", $Hopt, $maxQ, $Hmid_half_max);
  $cluster_max_agmr = $Hmid_half_max;
}
#print "ZZZZ\n";
#exit;

# construct graph with edges between vertices (accessions) for pairs with agmr < $cluster_max_agmr
# graph is constructed by adding edges and their endpoints; graph does not contain single unconnected vertices.
my $g = Graph::Undirected->new;
# keys: vertex pairs, values: agmrs


while (my ($e, $w) = each %edge_weight) {
 
  if ($w <= $cluster_max_agmr) {
    my ($id1, $id2) = split(" ", $e);
    $g->add_weighted_edge($id1, $id2, $w);
  }
}

my @ccs = $g->connected_components; # the connected components of graph are the clusters

my @output_lines = ();

my $count = 0;			# counts the number of accessions
my $count_cluster_accessions_out = 0;
for my $acc (@ccs) { # for each connected component (cluster of near-identical accessions)
  my %clusterids = map(($_ => 1), @$acc);
  my $cluster_noncluster_gap = least_noncluster_agmr(\%clusterids, \%id_closeidas);
  my $cc_size = scalar @$acc;
  $count += $cc_size;
  my $output_line_string = '';
  my @output_id_degree_pairs = ();
  my @sorted_cc = sort {$a cmp $b} @$acc; # sort the accession ids in the cluster
  my ($ccminw, $ccmaxw, $ccsumw, $nbad, $nbaddish) = (100000, -1, 0, 0, 0);
  while (my($i, $v) = each @sorted_cc) { # loop over every pair of ids in the cluster.
    my ($minw, $maxw) = (100000, -1);
    for (my $j=$i+1; $j<scalar @sorted_cc; $j++) {
      my $u = $sorted_cc[$j];
      if ($u eq $v) {
	warn "# $u $v  Why are they the same?\n";
	next;
      }
      my $edge_verts = ($v lt $u)? "$v $u" : "$u $v";
      my $weight = $edge_weight{$edge_verts} // -1;
      if ($weight == -1) {
	$nbad++;    # just counts pairs in cluster with agmr not found
      } else {
	$minw = $weight if($weight < $minw);
	$maxw = $weight if($weight > $maxw);
	$ccsumw += $weight;
	$nbaddish++ if($weight > $cluster_max_agmr)
      }
    }
    $ccminw = $minw if($minw < $ccminw);
    $ccmaxw = $maxw if($maxw > $ccmaxw);
    $output_line_string .= "$v ";
    my $degree = $g->degree($v);
    $output_line_string .= "$degree  ";
    push @output_id_degree_pairs, [$v, $degree];
  }
  @output_id_degree_pairs = sort {$b->[1] <=> $a->[1]} @output_id_degree_pairs; # sort by degree high to low
  my @id_degree_strs = map($_->[0] . ' ' . $_->[1], @output_id_degree_pairs);
  $output_line_string = join("  ", @id_degree_strs);
  my $n_cluster_edges = $cc_size*($cc_size-1)/2;
  my $ccavgw = $ccsumw/$n_cluster_edges;
  $output_line_string = sprintf("%4d  %7.5f %7.5f %7.5f %7.5f  %4d %4d  %s\n ",
				$cc_size, $ccminw, $ccavgw, $ccmaxw, $cluster_noncluster_gap,
				$nbaddish, $nbad, $output_line_string);
  if($cluster_noncluster_gap/$ccmaxw > 3){
    push @output_lines, $output_line_string;
    $count_cluster_accessions_out += $cc_size;
  }
}

open my $fhout, ">", "$output_cluster_filename" or die "Couldn't open $output_cluster_filename for writing.\n";
print $fhout "#  size d_min d_avg d_max d_cnc_min nbad1 nbad2 \n";
my @sorted_output_lines = sort { compare_str($a, $b) }  @output_lines;
print $fhout "# graph max edge length: $cluster_max_agmr. Found ", scalar @output_lines, " groups, with total of $count_cluster_accessions_out accessions.\n";
print $fhout join('', @sorted_output_lines);
close $fhout;

sub compare_str{ # sort by size of cluster, tiebreaker is avg intra-cluster distance.
  my $str1 = shift;
  my $str2 = shift;
  my @cols1 = split(" ", $str1);
  my @cols2 = split(" ", $str2);
  return ($cols1[0] != $cols2[0])? $cols1[0] <=> $cols2[0] : $cols1[2] cmp $cols2[2];
}

sub least_noncluster_agmr{	# call once for each cluster,
  # to get least agmr from cluster to outside cluster.
  my $clustids = shift;		# hash ref, keys ids in cluster.
  my $id1_id2as = shift; # hash ref; key ids; value: array ref of strings with id2, agmr12
  # my $cluster_max_agmr = shift;
  my $min_agmr_to_noncluster = 1;
  for my $id1 (keys %$clustids) { # loop over elements of cluster

    for my $s (@{$id1_id2as->{$id1}}) { # loop over accessions with smallish agmr to id1
      my ($id2, $d) = split(" ", $s);
      if (
	  # $d > $cluster_max_agmr  and  #
	  !exists $clustids->{$id2}
	 ) {
	$min_agmr_to_noncluster = $d if($d < $min_agmr_to_noncluster);
      }
    }
  }
  return $min_agmr_to_noncluster;
}
