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

# read duplicatesearch output and find clusters of accessions
# with near-identical genotypes
# specifically, create a graph, and for any pairs in duplicatesearch output with
# distance <= $link_max_distance, make an edge in the graph between the accessions of the pair.
# the clusters are then the connected components of the graph.
# i.e. two accession belong to the same cluster iff you can get from
# one to the other by traversing edges of the graph.
# one line of output for each cluster, with fields
# cluster size, min distance, max distance, 'nbad (should be 0)', the ids of accessions in cluster.
# e.g.:    3  0.0104 0.0142 0  TMEB778:250254008  TMEB778:250304613  TMEB778:250597946
# the duplicatesearch output should have all the distances for all the pairs in each cluster,
# and in that case the 4th field ('nbad') will be 0. If nbad > 0,
# consider rerunning this script with a smaller $max_distance,
# or rerun duplicatesearch with a larger value of 'max_estimated_distance' (-e option)

# usage example:
# distance_cluster.pl -in duplicatesearch.out  -out  distance_clusters  [-cluster_max 0.08]

# runs duplicatesearch, then distance_cluster, and outputs a file
# with the same format as duplicatesearch input, but now with just one
# line representing each cluster.

my $distances_filename = undef;
my $genotypes_filename = undef;
my $link_max_distance = 'auto'; # construct graph with edges between pairs of accessions iff their distance is <= this.
my $maxD = 1.0; # just ignore pairs separated by greater distance than this.
my $output_cluster_filename = "distance_cluster.out";
my $pow = 1;			# 'log';
my $minx = 0.001;
my $column = 6; # the distances to use for clustering are found in this column (unit-based)
my $in_out_factor = 2;
my $min_minextra_maxintra_ratio = -1; # can use this to only output clusters whose
# ratio ( min distance to extra-cluster accession / max intra-cluster distance )
# is at least this large.

GetOptions(
	   'distances_file=s' => \$distances_filename, # file with id1 id2 x xx distance_est distance (duplicatesearch output)
	   'genotypes_file|gt_file=s' => \$genotypes_filename,
	   'output_file=s' => \$output_cluster_filename,
	   'max_link_distance=f' => \$link_max_distance, # cluster using graph with edges for pairs with distance < this.
	   'pow=s' => \$pow,
	   'column=i' => \$column,
	   'min_ratio=f' => \$min_minextra_maxintra_ratio,
	  );
print STDERR "Clustering based on column: $column\n";
# $column--;
if (!defined $distances_filename) {
  print STDERR "Basic usage example: \n", "distance_cluster -in duplicatesearch.out  -out acluster.out \n";
  print STDERR "by default distance_cluster will attempt to automatically decide the max distance between duplicates.\n",
    " but you can specify it with the cluster option, e.g.  -cluster 0.06 \n";
  exit;
}


# store all distances in the distances file in hash %edge_weight
my %edge_weight = (); # keys are ordered pairs of accession ids representing graph edges, values are distances
my %id_closeidds = (); # keys are ids, values array ref of array refs of ids and distances of other accessions.
# all pairs found in duplicatesearch output are included. id1:[[$id2, $distance12], [$id3, $distance13], ...]
# if(0){
#   $column--;
# #my %id_distances = (); # keys are ids, values array ref of distances found between that id and others
# open my $fhin, "<", "$distances_filename" or die "Couldn't open $distances_filename for reading.\n";
# while (my $line = <$fhin>) {
#   next if($line =~ /^\s*#/);
#   my @cols = split(" ", $line);
#   my ($id1, $id2) =
#     # , $usable_chunks, $match_chunks, $est_distance) =
#     @cols[0..1];
#   my $distance = $cols[$column];
#   next if($distance > $maxD);
#   my $edge_verts = ($id1 lt $id2)? "$id1 $id2" : "$id2 $id1"; # order the pair of ids
#   if (!exists $id_closeidds{$id1}) {
#     $id_closeidds{$id1} = [[$id2, $distance]];
#   } else {
#     push @{$id_closeidds{$id1}}, [$id2, $distance];
#   }
#   # if(!exists $id_distances{$id1}){
#   #   $id_distances{$id1} = [$distance];
#   # }else{
#   #   push @{$id_distances{$id1}}, $distance;
#   # }
#   $edge_weight{$edge_verts} = $distance;
# }
# close $fhin;
# }else{
  my ($ew, $icid) = store_distances($distances_filename, $column);
  %edge_weight = %$ew;
  %id_closeidds = %$icid;
#}

# while(my($id1, $closeidds) = each %id_closeidds){
#   my @sort_closeidds = sort {$a->[1] <=> $b->[1]} @$closeidds;
#    print STDERR "#:  ", join(", ", map($_->[1], @sort_closeidds)), "\n";
#   $id_closeidds{$id1} = \@sort_closeidds;
 
# }


my @distances = values %edge_weight;
if ($link_max_distance eq 'auto') {
 
  my $cluster1d_obj = Cluster1d->new({label => '', xs => \@distances, pow => $pow, minx => $minx});
  print  "before two_cluster to choose cluster max distance\n";
  print  "#  pow: $pow  min: $minx \n";
  my ($Hopt, $maxQ, $Hmid_half_max) = $cluster1d_obj->two_cluster();
  # my ($Hoptx, $maxQx) = (0, 0); # $cluster1d_obj->two_cluster_x();
  $link_max_distance = $Hmid_half_max;
  print "# after two_cluster. Cluster max distance: $link_max_distance \n";
  #  print "# before k-means/kde clustering \n";
  #  my ($n_pts, $km_n_L, $km_n_R, $km_h_opt, $q, $kde_n_L, $kde_n_R, $kde_h_opt, $kde_q) = $cluster1d_obj->one_d_2cluster();
  #  printf( STDERR "# clustering %5d points;  k-means: %5d below  %8.6f and  %5d above; q: %6.4f.  kde: %5d below  %8.6f  and %5d above; kde_q: %6.4f   Hopt: %6.4f  maxQ: %6.4f.  Hmhmx: %6.4f \n",
  #       $n_pts, $km_n_L, $km_h_opt, $km_n_R, $q, $kde_n_L, $kde_h_opt, $kde_n_R, $kde_q, $Hopt, $maxQ, $Hmid_half_max);
  printf(STDERR  "Hopt: %6.4f  maxQ: %6.4f.  Hmhmx: %6.4f \n", $Hopt, $maxQ, $Hmid_half_max);
  
}
#print "ZZZZ\n";
#print STDERR "link max distance: $link_max_distance \n";
#exit;

# construct graph with edges between vertices (accessions) for pairs with distance < $link_max_distance
# graph is constructed by adding edges and their endpoints; graph does not contain single unconnected vertices.
my $g = Graph::Undirected->new;
# keys: vertex pairs, values: distances


while (my ($e, $w) = each %edge_weight) {
 
  if ($w <= $link_max_distance) {
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
  my ($cluster_noncluster_gap, $N_nearby_noncluster_pts, $N_cluster_pts_with_nearby_noncluster_pt) =
    least_noncluster_distance(\%clusterids, \%id_closeidds, $in_out_factor*$link_max_distance);
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
	$nbad++;    # just counts pairs in cluster with distance not found
      } else {
	$minw = $weight if($weight < $minw);
	$maxw = $weight if($weight > $maxw);
	$ccsumw += $weight;
	$nbaddish++ if($weight > $link_max_distance); # count intra-cluster pairs with distance > $link_max_distance
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
  $output_line_string = sprintf("%4d  %7.5f %7.5f %7.5f %7.5f  %4d %4d %4d %4d  %s\n ",
				$cc_size, $ccminw, $ccavgw, $ccmaxw, $cluster_noncluster_gap,
				$N_nearby_noncluster_pts, $N_cluster_pts_with_nearby_noncluster_pt,
				$nbaddish, $nbad, $output_line_string);
  if($cluster_noncluster_gap/$ccmaxw >= $min_minextra_maxintra_ratio){
    push @output_lines, $output_line_string;
    $count_cluster_accessions_out += $cc_size;
  }
}

open my $fhout, ">", "$output_cluster_filename" or die "Couldn't open $output_cluster_filename for writing.\n";
print $fhout "# graph max edge length: $link_max_distance. Found ", scalar @output_lines, " groups, with total of $count_cluster_accessions_out accessions.\n";
print $fhout "# col 1: cluster size.\n";
print $fhout "# col 2, 3, 4: min, avg. and max. intra-cluster distances.\n";
print $fhout "# col 5: min distance between cluster pts and non-cluster pts.\n";
print $fhout "# col 6: number of 'nearby' non-cluster pts.\n";
print $fhout "# col 7: number of cluster pts with 'nearby' non-cluster pts.\n";
print $fhout "# 'nearby' defined as within $in_out_factor * link_max_distance = ", $in_out_factor*$link_max_distance, ".\n";
print $fhout "# col 8: number of intra-cluster distance > link_max_distance ($link_max_distance).\n";
print $fhout "# col 9: number of intra-cluster pairs with distance not present in input file.\n";
print $fhout "# then clusters and for each the number of cluster members withing link_max_distance.\n";
my @sorted_output_lines = sort { compare_str($a, $b) }  @output_lines;
print $fhout join('', @sorted_output_lines);
close $fhout;

sub compare_str{ # sort by cluster size; 1st tiebreaker: avg intra-cluster distance; 2nd tiebreaker: first accession id. 
  my $str1 = shift;
  my $str2 = shift;
  my @cols1 = split(" ", $str1);
  my @cols2 = split(" ", $str2);
  return ($cols1[0] != $cols2[0])? $cols1[0] <=> $cols2[0] :   # cluster size.
    ($cols1[2] <=> $cols2[2])? $cols1[2] <=> $cols2[2] :  # avg. intra-cluster distance.
    $cols1[7] cmp $cols2[7]; #  id of first accession in cluster.
}

sub least_noncluster_distance{	# call once for each cluster,
  # to get least distance from cluster to outside cluster.

  my $clustids = shift;		# hash ref, keys ids in cluster.
  my $id1_id2ds = shift; # hash ref; key ids; value: array ref of strings with id2, distance12
  my $near_noncluster_distance = shift; # set to some multiple of $link_max_distance
  my $short_cluster_noncluster_edge_count = 0; # count number of cluster-noncluster edges < $near_noncluster_distance
  my %near_noncluster_points = ();
  my %cluster_pts_near_noncluster_pt = ();
  my $min_distance_to_noncluster = 1;
  for my $id1 (keys %$clustids) { # loop over elements of cluster
    for my $s (@{$id1_id2ds->{$id1}}) { # loop over accessions with smallish distance to id1
      my ($id2, $d) = @$s;
      if ( !exists $clustids->{$id2} ) { # $id2 is not in the cluster
	if($d < $min_distance_to_noncluster){
	  $min_distance_to_noncluster = $d;
	}
	#	else{ # @$ is sorted by $d, small to large, so if this $d isn't small enough, we're done with @$s
	if($d < $near_noncluster_distance){
	  $short_cluster_noncluster_edge_count++;
	  $near_noncluster_points{$id2}++;
	  $cluster_pts_near_noncluster_pt{$id1}++;
	}else{
	  last; #
	}
#	}
      }
    }
  }
  my $nearby_noncluster_points_count = scalar keys %near_noncluster_points;
  my $cluster_pts_near_noncluster_pt_count = scalar keys %cluster_pts_near_noncluster_pt;
  return ($min_distance_to_noncluster, $nearby_noncluster_points_count, $cluster_pts_near_noncluster_pt_count);
}

sub store_distances{
  my $distances_filename = shift;
  my $column = shift; # (unit based) column number in which to find the distances.
  print STDERR "column: $column\n";
  $column--;
my %edge_weight = (); # keys are ordered pairs of accession ids representing graph edges, values are distances
my %id_closeidds = (); # keys are ids, values array ref of array refs of ids and distances of other accessions.
# all pairs found in duplicatesearch output are included. id1:[[$id2, $distance12], [$id3, $distance13], ...]
#my %id_distances = (); # keys are ids, values array ref of distances found between that id and others
open my $fhin, "<", "$distances_filename" or die "Couldn't open $distances_filename for reading.\n";
while (my $line = <$fhin>) {
  next if($line =~ /^\s*#/);
  my @cols = split(" ", $line);
  my ($id1, $id2) =
    # , $usable_chunks, $match_chunks, $est_distance) =
    @cols[0..1];
  my $distance = $cols[$column];
  next if($distance > $maxD);
  my $edge_verts = ($id1 lt $id2)? "$id1 $id2" : "$id2 $id1"; # order the pair of ids
  if (!exists $id_closeidds{$id1}) {
    $id_closeidds{$id1} = [[$id2, $distance]];
  } else {
    push @{$id_closeidds{$id1}}, [$id2, $distance];
  }
  # if(!exists $id_distances{$id1}){
  #   $id_distances{$id1} = [$distance];
  # }else{
  #   push @{$id_distances{$id1}}, $distance;
  # }
  $edge_weight{$edge_verts} = $distance;
}
  close $fhin;
  return (\%edge_weight, \%id_closeidds);
}
