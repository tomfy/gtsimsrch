#!/usr/bin/perl -w
use strict;
use Graph::Undirected;
use Getopt::Long;
use File::Basename 'dirname';
use List::Util qw 'min max sum';

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

# read duplicate_search output and find clusters of accessions
# with near-identical genotypes
# specifically, create a graph, and for any pairs in duplicate_search output with
# distance <= $link_max_distance, make an edge in the graph between the accessions of the pair.
# the clusters are then the connected components of the graph.
# i.e. two accession belong to the same cluster iff you can get from
# one to the other by traversing edges of the graph.
# one line of output for each cluster, with fields
# cluster size, min distance, max distance, 'missing_distance_count (should be 0)', the ids of accessions in cluster.
# e.g.:    3  0.0104 0.0142 0  TMEB778:250254008  TMEB778:250304613  TMEB778:250597946
# the duplicate_search output should have all the distances for all the pairs in each cluster,
# and in that case the 4th field ('missing_distance_count') will be 0. If missing_distance_count > 0,
# consider rerunning this script with a smaller $max_distance,
# or rerun duplicate_search with a larger value of 'max_estimated_distance' (-e option)

# usage example:
# clusterer.pl -in duplicate_search.out  -out  distance_clusters  [-cluster_distance 0.08]

# runs duplicate_search, then distance_cluster, and outputs a file
# with the same format as duplicate_search input, but now with just one
# line representing each cluster.

{				# beginning of 'main'
  my $command_str = "# command:  $0  " . join(" ", @ARGV);
  
  my $distances_filename = undef;
  my $link_max_distance = 'auto'; # construct graph with edges between pairs of accessions iff their distance is <= this.
  my $maxD = 1; # just ignore pairs separated by greater distance than this.
  my $output_cluster_filename = "distance_cluster.out";
  my $pow = 1; # cluster transformed values tx_i = pow(x_i, $pow), or if $pow is 'log' tx_i = log(x_i)
# column numbers are unit-based:
  my $id1_column = 1;
  my $mdc1_column = 2; # missing data count, accession 1
  my $id2_column = 3;
  my $mdc2_column = 4; # missing data count, accession 2
  my $d_column = 5; # the distances to use for clustering are found in this column (unit-based)
  my $in_out_factor = 1.2;
  my $f = 0.2; # fraction of way auto link_max_distance is across the Q > 0.5*Qmax range.
  my $full_output = 0; # if '-nofull' will not output the degrees, etc for each cluster member, etc.
  my $prune_factor = 0.2; # if cluster member has edges to < this fraction of other cluster members, delete it.

  my $min_minextra_maxintra_ratio = -1; # can use this to only output clusters whose
  # ratio ( min distance to extra-cluster accession / max intra-cluster distance )
  # is at least this large.

  GetOptions(
	     'distances_file|input=s' => \$distances_filename, # file with id1 id2 distance  ... (duplicate_search output)
	     'cluster_distance|link_distance|dlink=s' => \$link_max_distance, # cluster using graph with edges for pairs with distance < this.
	     'output_file=s' => \$output_cluster_filename,

	     'maxd|dmax=f' => \$maxD,
	     'pow=s' => \$pow,
	     'id1col=i' => \$id1_column,
	     'id2col=i' => \$id2_column,
	     'dcolumn=i' => \$d_column,
	     'min_ratio=f' => \$min_minextra_maxintra_ratio,
	     'prune_factor=f' => \$prune_factor,
	     'full_output!' => \$full_output,
	    );

  if (!defined $distances_filename) {
    print STDERR "Input file must be specified.\n";
    print STDERR "Basic usage example: \n", "clusterer -in duplicate_search.out  -out cluster.out \n";
    print STDERR "by default distance_cluster will attempt to automatically decide the max distance between duplicates.\n",
      " but you can specify it with the cluster_distance option, e.g.  -cluster_distance 0.06 \n";
    exit;
  }
  
  print "# Getting accessions ids from columns: $id1_column, $id2_column \n";
  print "# Clustering based on distances in column: $d_column\n";
  print STDERR "Pair distances filename: $distances_filename \n";

  #######################################################################################################
  # store all distances in the distances file in hash refs $edge_weight and id_closeidds
  # $edge_weight : keys are ordered pairs of accession ids representing graph edges, values are distances
  # $id_closeidds : keys are ids, values array ref of array refs of ids and distances of other accessions.

  #print STDERR "$distances_filename $id1_column $id2_column $d_column $maxD\n";
  my ($edge_weight, $id_closeidds, $id_mdcount) = store_distances($distances_filename,
								  $id1_column, $mdc1_column, $id2_column, $mdc2_column,
								  $d_column, $maxD);
  #######################################################################################################
  # get array of edges ordered by weight (distance), small to large.
  my @sorted_edges = sort {$edge_weight->{$a} <=> $edge_weight->{$b} } keys %$edge_weight;
  print "# Number of edges stored ", scalar @sorted_edges, "\n";

  #######################################################################################################
  # if $link_max_distance not specified, attempt to find a reasonable value by looking at the distances
  # found in the distances input file.
#  my ($Hopt, $Qmax, $Ledge, $Redge, $H, $Q);
  if ($link_max_distance eq 'auto') {
  my  ($Hopt, $Qmax, $Ledge, $Redge, $H, $Q) = auto_max_link_distance($edge_weight, $pow, $f);
    $link_max_distance = $H;
    print "#  $H  $Q \n";
    #  print "## $Hoptxx  $Qmaxxx  $Hxx  $Qxx \n";
  }
  print "# Max link distance: $link_max_distance\n";
  #######################################################################################################
  

  #######################################################################################################
  # construct graph with edges between vertices (accessions) for pairs with distance < $link_max_distance
  # construct graph by adding edges & their endpoints; graph does not contain single unconnected vertices.

  my $the_graph = construct_graph($edge_weight, $link_max_distance, \@sorted_edges);
  # the connected components of graph are the clusters
  my @clusters = $the_graph->connected_components; # array of array refs of ids
  my @vertices = $the_graph->vertices();
  my $n_nodes = scalar @vertices;
  print "# Graph created with ", scalar @clusters, " clusters and $n_nodes cluster members.\n";
  ########################################################################################################

  prune_graph($the_graph, $prune_factor);
  @clusters = $the_graph->connected_components; # array of array refs of ids
  @vertices = $the_graph->vertices();
  $n_nodes = scalar @vertices;
  ########################################################################################################
  # loop over clusters, evaluating quality of each.
  
  my @output_lines = ();
  my $count_cluster_accessions = 0; # counts the number of accessions in clusters (should be equal to $n_nodes)
  my $count_cluster_accessions_out = 0;
  my $count_rejected_clusters = 0;
  while ( my($i, $accs) = each @clusters) { # for each connected component (cluster of near-identical accessions)
    my $cluster_size = scalar @$accs;
    $count_cluster_accessions += $cluster_size;
    my %thisclusterids = map(($_ => 1), @$accs);

    my ($cluster_min_d, $cluster_avg_d, $cluster_max_d, $intracluster_far_pair_count, $missing_distance_count, $out_iddegds, $min_degree) =
      cluster_quality1($accs, $edge_weight, $the_graph, $link_max_distance);

    my ($cluster_noncluster_gap, $N_nearby_noncluster_pts, $N_cluster_pts_with_nearby_noncluster_pt) =
      cluster_quality2($the_graph, \%thisclusterids, $id_closeidds, $link_max_distance, $in_out_factor*$link_max_distance);

    my $cluster_q_category = 0;
    $cluster_q_category++ if($intracluster_far_pair_count > 0);
    $cluster_q_category++ if($cluster_max_d > 1.2*$link_max_distance);
    $cluster_q_category += min($N_nearby_noncluster_pts, $N_cluster_pts_with_nearby_noncluster_pt, 2);
    $cluster_q_category++ if($cluster_noncluster_gap < 1.2*$link_max_distance);
    
    my $output_line_string = '';
    if ($full_output) { # for each duplicate group member output  id, degree, rms dist from duplicates
      # sort by missing data, low to high.
      my @sorted_iddegds = sort {$id_mdcount->{$a->[0]} <=> $id_mdcount->{$b->[0]}} @$out_iddegds;
      
      $output_line_string = join("  ", map(  sprintf("%s %1d %7.5f  ", $_->[0], $_->[1], $_->[3]) ,@sorted_iddegds) );
    } else {
      #my @ids = map( sprintf("%s  ", $_->[0]), @$out_iddegds);
      my @ids = map($_->[0], @$out_iddegds);
      my @sorted_ids = sort {$id_mdcount->{$a} <=> $id_mdcount->{$b}} @ids;
      print join(" ", map(sprintf(" %s %ld ", $_, $id_mdcount->{$_}), @sorted_ids)), "\n";
      $output_line_string = join("  ", @sorted_ids);
    }
    $output_line_string = sprintf("%4d  %7.5f %7.5f %7.5f %7.5f  %4d %4d %4d %4d   %s\n ",
				  $cluster_size, $cluster_min_d, $cluster_avg_d, $cluster_max_d, $cluster_noncluster_gap,
				  $N_nearby_noncluster_pts, $N_cluster_pts_with_nearby_noncluster_pt,
				  $intracluster_far_pair_count+$missing_distance_count, $min_degree,
				  #  $d, $e,
				  #	  $cluster_q_category,
				  $output_line_string);
    if ($cluster_max_d == 0  or  $cluster_noncluster_gap/$cluster_max_d >= $min_minextra_maxintra_ratio) {
      push @output_lines, $output_line_string;
      $count_cluster_accessions_out += $cluster_size;
    } else {
      $count_rejected_clusters++;
    }
  }				# end loop over clusters
  #######################################################################################################

  ########################################################################################################
  # output
  
  my $N_clusters_out = scalar @output_lines;
  die if($N_clusters_out + $count_rejected_clusters != scalar @clusters);
  print "# Found $N_clusters_out acceptable clusters containing $count_cluster_accessions_out accessions.\n";
  print "# Removing duplicates will reduce the number of accessions by ", $count_cluster_accessions_out - $N_clusters_out , "\n"; 
  open my $fhout, ">", "$output_cluster_filename" or die "Couldn't open $output_cluster_filename for writing.\n";
  print $fhout $command_str, "\n";
  print $fhout "# Graph max edge length: $link_max_distance. \n";
  print $fhout "# Found $N_clusters_out acceptable clusters, containing $count_cluster_accessions_out accessions.\n";
  print $fhout "# Removing duplicates will reduce the number of accessions by ", $count_cluster_accessions_out - $N_clusters_out , "\n"; 
  print $fhout "# col 1: cluster size.\n";
  print $fhout "# col 2, 3, 4: min, avg. and max. intra-cluster distances.\n";
  print $fhout "# col 5: min distance between cluster pts and non-cluster pts.\n";
  print $fhout "# col 6: number of 'nearby' non-cluster pts.\n";
  print $fhout "# col 7: number of cluster pts with 'nearby' non-cluster pts.\n";
  print $fhout "# 'nearby' defined as within $in_out_factor * link_max_distance = ", $in_out_factor*$link_max_distance, ".\n";
  print $fhout "# col 8: number of intra-cluster distances > link_max_distance ($link_max_distance).\n";
  print $fhout "# col 9: number of intra-cluster pairs with distance not present in input file.\n";
  #  print $fhout "# col 10: max distance between cluster pt. and cluster consensus.\n";
  #  print $fhout "# col 11: number of cluster pts. further than $d_to_consensus_factor * link_max_distance from cluster consensus.\n";
  print $fhout "# then ids of cluster members, sorted with least missing data first.\n";
  print $fhout "# if -full option invoked, then each id is followed by:\n";
  print $fhout "#   degree of the accession and its rms distance from others in cluster.\n";
  print $fhout "#   and for each the number of other cluster members within link_max_distance.\n";
  my @sorted_output_lines = sort { compare_str($a, $b) }  @output_lines;
  print $fhout join('', @sorted_output_lines);
  close $fhout;

}				# end of 'main'

##########################################################################################################

sub compare_str{    # sort by cluster size;
  # 1st tiebreaker: avg intra-cluster distance;
  # 2nd tiebreaker: max intra-cluster distance,
  # 3rd tiebreaker min. cluster-noncluster distance,
  # 4th tiebreaker: first accession id.
  my $str1 = shift;
  my $str2 = shift;
  my @cols1 = split(" ", $str1);
  my @cols2 = split(" ", $str2);
  return ($cols1[0] != $cols2[0])? $cols1[0] <=> $cols2[0] : # cluster size, small to large.
    ($cols1[2] <=> $cols2[2])? $cols1[2] <=> $cols2[2] : # avg. intra-cluster distance, small to large.
    ($cols1[3] <=> $cols2[3])? $cols1[3] <=> $cols2[3] : # max. intra-cluster distance, small to large.
    ($cols2[4] <=> $cols1[4])? $cols2[4] <=> $cols1[4] : # min. cluster-noncluster distance, large to small.
    (defined $cols1[9])? ($cols1[9] cmp $cols2[9]) : 0;	#  id of first accession in cluster.
}

sub store_distances{
  my $distances_filename = shift;
  my $id1_col = shift;
  my $mdc1_col = shift;
  my $id2_col = shift;
  my $mdc2_col = shift;
  my $d_column = shift; # (unit based) column number in which to find the distances (agmrs).
  my $maxD = shift;
  $id1_col--; $mdc1_col--; $id2_col--; $mdc2_col--; $d_column--; # make them zero-based
  my $n_def = 0;
  my $n_undef = 0;
  my %edge_weight = (); # keys are ordered pairs of accession ids representing graph edges, values are distances
  my %id_closeidds = (); # keys are ids, values array ref of array refs of ids and distances of other accessions.
  my %id_mdc = (); # keys: ids; values: missing data counts.
  # all pairs found in duplicate_search output are included. id1:[[$id2, $distance12], [$id3, $distance13], ...]
  open my $fhin, "<", "$distances_filename" or die "Couldn't open $distances_filename for reading.\n";
  while (my $line = <$fhin>) {
    next if($line =~ /^\s*#/);
    my @cols = split(" ", $line);
    my ($id1, $mdc1, $id2, $mdc2) = @cols[$id1_col, $mdc1_col, $id2_col, $mdc2_col];
    my $distance = $cols[$d_column];
    $id_mdc{$id1} = $mdc1;
    $id_mdc{$id2} = $mdc2;
    #   print STDERR "$id1 $id2  $distance  $maxD\n";
    next if($distance > $maxD);
    my $edge_verts = ($id1 lt $id2)? "$id1 $id2" : "$id2 $id1"; # order the pair of ids
    if(!defined $edge_verts){
      print "id1: ", (defined $id1)? $id1 : 'undef', "\n";
      print "id2: ", (defined $id2)? $id2 : 'undef', "\n";
      $n_undef++;
    }else{
      $n_def++;
    }
    if (!exists $id_closeidds{$id1}) {
      $id_closeidds{$id1} = [[$id2, $distance]];
    } else {
      push @{$id_closeidds{$id1}}, [$id2, $distance];
    }
    $edge_weight{$edge_verts} = $distance;
  }
  close $fhin;
  print "n def, undef: $n_def $n_undef\n";
  return (\%edge_weight, \%id_closeidds, \%id_mdc);
}

sub auto_max_link_distance{
  my $edge_weight = shift;
  my $pow = shift;
  my $f = shift;
  my @distances = values %$edge_weight;
  # print STDERR "## number of distances: ", scalar @distances, "\n";

  my $cluster1d_obj = Cluster1d->new({label => '', xs => \@distances, pow => $pow});
#  print  "#  before two_cluster to choose cluster max distance\n";
#  print  "#  pow: $pow \n";
  my ($Hopt, $maxQ, $Hopt_avg, $maxQ_avg, $Ledge, $Redge) = $cluster1d_obj->find_cluster_at_left($f);
#  my $Hmid_half_max = 0.5*($Ledge+$Redge);
#  my $H33 = (0.67*$Ledge + 0.33*$Redge);
  #  printf(STDERR  "# Hopt: %6.4f  maxQ: %6.4f. Hoptavg maxQavg: %6.4f %6.4f   Hmhmx: %6.4f H33: %6.4f \n", $Hopt, $maxQ, $Hopt_avg, $maxQ_avg, $Hmid_half_max, $H33);
  print STDERR "$Hopt $maxQ  $Ledge $Redge  $Hopt_avg  $maxQ_avg \n";
  return ($Hopt, $maxQ, $Ledge, $Redge, $Hopt_avg, $maxQ_avg);
}

sub construct_graph{
  my $edge_weight = shift; # hash ref; keys: (ordered) id pairs; values edge weight (distance)
  my $max_link_distance = shift;
  my $sorted_edges = shift; # sorted by weight (distance) small to large
  print "number of keys of edge_weight: ", scalar keys %$edge_weight, "\n";
  print "number of sorted edges: ", scalar @$sorted_edges, "\n";
  my $the_graph = Graph::Undirected->new;
    for my $e (@$sorted_edges){
    my $w = $edge_weight->{$e};
    if ($w <= $max_link_distance) {
      my ($id1, $id2) = split(" ", $e);
      $the_graph->add_weighted_edge($id1, $id2, $w);
    } else {
      last;
    }
  }
  return $the_graph;
}

sub prune_graph{
  # for each cluster, remove nodes (and assoc. edges)
  # if degree too small compared to cluster size
  my $the_graph = shift;
  my $prune_factor = shift;
  my $vert_deg = {};
  for my $v ($the_graph->vertices()){
    $vert_deg->{$v} = $the_graph->degree($v);
  }
  my @clusters = $the_graph->connected_components;
  for my $clstr (@clusters) {
    prune_cluster($the_graph, $clstr, $vert_deg, $prune_factor);
  }
}

sub prune_cluster{
  my $graph = shift;
  my $cluster = shift;		# array ref of
  my $vertex_degree = shift;
  my $prune_factor = shift;
  my $cluster_size = scalar @$cluster;
  my $min_degree = $prune_factor*($cluster_size - 1);
  my $prune_count = 0;
  for (my $i=0; $i < scalar @$cluster; $i++) {
    my $u = $cluster->[$i];
    my $u_dubious = ($vertex_degree->{$u} < $min_degree);
    if ($u_dubious) {
      $graph->delete_vertex($u); # also deletes associated edges
      $prune_count++;
    }
  }
}

sub cluster_quality1{
  my $accs = shift;		# array ref of ids of cluster members
  my $edge_weight = shift; # hashref; keys: order id pairs, values: weights (distances)
  my $graph = shift;
  my $link_max_distance = shift;
  # my $id_degree = shift;
  my $output_line_string = '';
  my @output_id_degree_maxds = ();
  my @sorted_clusterids = sort {$a cmp $b} @$accs; # sort the accession ids in the cluster
  my $cluster_size = scalar @sorted_clusterids;
  my ($cluster_min_d, $cluster_max_d, $cluster_sum_d, $missing_distance_count, $intracluster_far_pair_count) = (10000, -1, 0, 0, 0);
  my $sum_of_degrees = 0;
  my $cluster_edge_count = 0;
  my %id_maxd = map {$_ => -1} @sorted_clusterids;
  my %id_sumdsq = map  {$_ => 0} @sorted_clusterids;

  while (my($i, $v) = each @sorted_clusterids) { # loop over every pair of ids in the cluster.
    my $v_degree = $graph->degree($v); # $id_degree->{$v}; # degree of $v (=number of other nodes at dist <= $link_max_distance
    $sum_of_degrees += $v_degree;    # graph->degree($v);
    $output_line_string .= "$v  $v_degree  ";
    #  $id_degree{$v} = $degree;
    my $present_distance_count = 0;
    my ($minw, $maxw) = (10000, -1);

    for (my $j=$i+1; $j<scalar @sorted_clusterids; $j++) {
      my $u = $sorted_clusterids[$j];
      if ($u eq $v) {
	warn "# $u $v  Why are they the same?\n";
	next;
      }
      my $edge_verts = ($v lt $u)? "$v $u" : "$u $v";
      my $weight = $edge_weight->{$edge_verts} // -1;
      if ($weight == -1) {
	$missing_distance_count++; # just counts pairs in cluster with distance not found
      } else {
	$present_distance_count++;
	$minw = $weight if($weight < $minw);
	$maxw = $weight if($weight > $maxw);
	$id_maxd{$v} = max($id_maxd{$v}, $weight);
	$id_maxd{$u} = max($id_maxd{$u}, $weight);
	$id_sumdsq{$v} += $weight*$weight;
	$id_sumdsq{$u} += $weight*$weight;

	$cluster_sum_d += $weight;
	$cluster_edge_count++;
	$intracluster_far_pair_count++ if($weight > $link_max_distance); # count intra-cluster pairs with distance > $link_max_distance
      }
    }			 # loop over $u (only those with indices > $i)
    $cluster_min_d = $minw if($minw < $cluster_min_d);
    $cluster_max_d = $maxw if($maxw > $cluster_max_d);

  } # loop over $v
  my $min_degree = scalar @sorted_clusterids;
  for my $id (@sorted_clusterids) {
    my $degree = $graph->degree($id); # $id_degree->{$id};
    $min_degree = $degree if($degree < $min_degree);
    my $rms_dist = sqrt($id_sumdsq{$id}/($cluster_size-1));
    push @output_id_degree_maxds, [$id, $degree, $id_maxd{$id}, $rms_dist];
  }  
  @output_id_degree_maxds = sort {$b->[1] <=> $a->[1]} @output_id_degree_maxds;
  return ($cluster_min_d, $cluster_sum_d/$cluster_edge_count, $cluster_max_d, $intracluster_far_pair_count, $missing_distance_count, \@output_id_degree_maxds, $min_degree)
}

sub cluster_quality2{
  my $graph = shift;
  my $thisclustids = shift;	# hash ref, keys ids in this cluster.
  my $id1_id2ds = shift; # hash ref; key ids; value: array ref of strings with id2, distance12
  my $max_link_distance = shift;
  my $near_noncluster_distance = shift; # set to some multiple of $link_max_distance

  my $min_distance_to_noncluster = 1;
  my $short_cluster_noncluster_edge_count = 0;
  my %near_noncluster_points = ();
  my %cluster_pts_near_noncluster_pt = ();
  my $cluster_size =  scalar keys %$thisclustids;
  my $sum_of_degrees = 0;
  my ($min_intracluster_d, $max_intracluster_d, $avg_intracluster_d) = (1000, -1, 0);
  my $cluster_edge_count = 0;
  for my $id1 (keys %$thisclustids) { # loop over elements of cluster
    my $degree = $graph->degree($id1);
    $sum_of_degrees += $degree;
    for my $s (@{$id1_id2ds->{$id1}}) { # loop over accessions with smallish distance to id1
      my ($id2, $d) = @$s;
      if ( exists $thisclustids->{$id2} ) { # $id2 is in the cluster
	$min_intracluster_d = min($d, $min_intracluster_d);
	$max_intracluster_d = max($d, $max_intracluster_d);
	$avg_intracluster_d += $d;
	$cluster_edge_count++;
      } else {			# $id2 is not in this cluster
	if ($d < $min_distance_to_noncluster) {
	  $min_distance_to_noncluster = $d;
	}
	if ($d < $near_noncluster_distance) {
	  $short_cluster_noncluster_edge_count++;
	  $near_noncluster_points{$id2}++;
	  $cluster_pts_near_noncluster_pt{$id1}++;
	} else {
	  last;			#
	}
      }			     # end id2 not in clustr
    }			     # end loop over nearish accessions to id1
  }
  $avg_intracluster_d /= ($cluster_size*($cluster_size-1));
  my $xxx = ($cluster_size*($cluster_size-1) - $sum_of_degrees)/2;
  my $nearby_noncluster_points_count = scalar keys %near_noncluster_points;
  my $cluster_pts_near_noncluster_pt_count = scalar keys %cluster_pts_near_noncluster_pt;
  return ($min_distance_to_noncluster, $nearby_noncluster_points_count, $cluster_pts_near_noncluster_pt_count)
}

################################################################################################################
################################################################################################################

# sub store_dosages{
#   my $dosages_filename = shift;
#   my $id_closeidds = shift; 
#   # my $clusterids = shift; # hash ref; keys are ids; values array refs of genotypes (0, 1, 2, ...)
#   # store individual lines of dosage file in hash
#   open my $fh_dosages, "<", "$dosages_filename";
#   #open my $fhout, ">", "$dosages_output_filename" or die "Couldn't open $dosages_output_filename for writing.\n";

#   # store ids and genotypes of clusters (size >= 2).
#   my %id_gts  = (); # key: ids; value: array ref of dosages (0,1,2, $missing_data_char) 
#   while (my $first_line = <$fh_dosages>) { 
#     next if($first_line =~ /^\s*#/); # skip comments
#     last;
#   }
#   while (my $line = <$fh_dosages>) {
#     next if($line =~ /^\s*#/);	# skip comments
#     my @cols = split(" ", $line);
#     my $id = shift @cols;
#     if (exists $id_closeidds->{$id}) {
#       #   if (exists $clusterids->{$id}) { # store if it is in a cluster
#       $id_gts{$id} = \@cols;
#       #   }
#     }
#   }
#   close($fh_dosages);
#   print STDERR "# Stored ", scalar keys %id_gts, " accessions and their genotypes.\n";
#   return \%id_gts;
# }

# sub cluster_consensus_dosages{ # get a consensus genotypes set for one cluster.
#   my $cluster_ids = shift; # array ref holding the ids of the accessions in the cluster.
#   my $id_dosages = shift; # hash ref. key: id, value: array ref of dosages for this id.
#   my $first_id = $cluster_ids->[0];
#   my $first_dosages = $id_dosages->{$first_id};
#   my $n_markers = scalar @$first_dosages;
#   my @consensus_dosages = (0) x $n_markers;

#   my $cluster_size = scalar @$cluster_ids;
#   if ($cluster_size > 2) {
#     my @marker_votes = ();
#     for (1..$n_markers) {
#       push @marker_votes, [0, 0, 0, 0]; # counts of dosages 0, 1, 2, and missing data
#     }
#     for my $an_id (@$cluster_ids) { # collect the genotype distributions for all markers
#       my $dosages = $id_dosages->{$an_id}; # array ref of dosages for $an_id
#       while (my($i, $d) = each @$dosages) {
# 	$d = 3 if($d eq $missing_data_char);
# 	$marker_votes[$i]->[$d]++;
#       }
#     }
#     while (my($i, $mv) = each @marker_votes) { # 
#       my $e = $missing_data_char;
#       #   print STDERR "    $i  ", join(" ", @$mv), "\n";
#       for my $j (0..2) {
# 	my $vote = $mv->[$j]; # the number of votes for dosage = $j for marker $i
# 	if (3*$vote >= 2*$cluster_size) { # must have at least 2/3 of one dosage, or call it missing data. 
# 	  $e = $j;
# 	  last;
# 	}
#       }
#       $consensus_dosages[$i] = $e;
#     }
#     #  if($cluster_size > 20){
#     # while(my($i, $mv) = each @marker_votes){
#     #   print STDERR "XXX: $cluster_size   $i  ", join("  ", @$mv), "  ", $consensus_dosages[$i], "\n";
#     # }
#     # }
#     #}
#   } else {			# cluster size == 2
#     while (my($i, $d1) = each @$first_dosages) {
#       my $second_id = $cluster_ids->[1];
#       my $d2 = $id_dosages->{$second_id}->[$i];
#       if ($d1 eq $d2) {
# 	$consensus_dosages[$i] = $d1;
#       } elsif ($d1 eq $missing_data_char) { # d1 is missing, use d2 (which might also be missing)
# 	$consensus_dosages[$i] = $d2;
#       } elsif ($d2 eq $missing_data_char) { # d2 is missing, use d1 (which is not missing)
# 	$consensus_dosages[$i] = $d1;
#       } else {	 # neither is missing, and they disagree - use average
# 	$consensus_dosages[$i] = 0.5*($d1 + $d2);
#       }
#     }
#   }
 
#   return \@consensus_dosages;
# }

# sub distance{
#   my $dosages1 = shift;
#   my $dosages2 = shift;
#   my ($numer, $denom) = (0, 0); 
#   while (my($i, $dsg1) = each @$dosages1) {
#     next if($dsg1 eq $missing_data_char);
#     my $dsg2 = $dosages2->[$i];
#     next if($dsg2 eq $missing_data_char);
#     $numer += abs($dsg1 - $dsg2);
#     $denom++;
#   }
#   # print STDERR "## $numer  $denom  ", $numer/$denom, "\n";
#   return ($denom > 0)? $numer/$denom : $bad_distance;
# }

# sub cluster_quality_consensus{
#   my $consensus_dosages = shift;
#   my $id_gts = shift;
#   my $thisclustids = shift;	# hash ref, keys ids in this cluster.
#   my $id1_id2ds = shift; # hash ref; key ids; value: array ref of strings with id2, distance12
#   my $max_link_distance = shift;
#   my $near_noncluster_distance = shift; # set to some multiple of $link_max_distance

#   my @cluster_distances_from_consensus = ();
#   my $max_d_to_consensus = -1;
#   my $far_from_consensus_count = 0;
#   my $min_distance_to_noncluster = 1;
#   my $short_cluster_noncluster_edge_count = 0;
#   my %near_noncluster_points = ();
#   my %cluster_pts_near_noncluster_pt = ();
#   my $cluster_size =  scalar keys %$thisclustids;
#   # print STDERR "# Cluster size: $cluster_size  ##  $near_noncluster_distance \n";
#   my %clusterpt_dist2consensus = ();
#   for my $id1 (keys %$thisclustids) { # loop over elements of cluster
#     my $d_to_consensus = distance($consensus_dosages, $id_gts->{$id1});
#     $max_d_to_consensus = max($d_to_consensus, $max_d_to_consensus);
#     if ($d_to_consensus > $d_to_consensus_factor*$max_link_distance ) {
#       $far_from_consensus_count++;
#     }
#     $clusterpt_dist2consensus{$id1} = $d_to_consensus;
#     for my $s (@{$id1_id2ds->{$id1}}) { # loop over accessions with smallish distance to id1
#       my ($id2, $d) = @$s;
#       #  print STDERR "#    id2: d(id1, id2):  $id2  $d \n";
#       if ( !exists $thisclustids->{$id2} ) { # $id2 is not in this cluster
# 	if ($d < $min_distance_to_noncluster) {
# 	  $min_distance_to_noncluster = $d;
# 	}
# 	#	else{ # @$ is sorted by $d, small to large, so if this $d isn't small enough, we're done with @$s
# 	if ($d < $near_noncluster_distance) {
# 	  $short_cluster_noncluster_edge_count++;
# 	  $near_noncluster_points{$id2}++;
# 	  $cluster_pts_near_noncluster_pt{$id1}++;
# 	} else {
# 	  last;			#
# 	}
#       }			     # end id2 not in clustr
#     }			     # end loop over nearish accessions to id1
#   }
#   #print STDERR "# $max_d_to_consensus  $far_from_consensus_count \n";
#   my $nearby_noncluster_points_count = scalar keys %near_noncluster_points;
#   my $cluster_pts_near_noncluster_pt_count = scalar keys %cluster_pts_near_noncluster_pt;
#   #  printf (STDERR "# %8.6f %5d %5d   %8.6f %5d\n", $min_distance_to_noncluster, $nearby_noncluster_points_count,
#   #	  $cluster_pts_near_noncluster_pt_count, $max_d_to_consensus, $far_from_consensus_count);
#   return ($min_distance_to_noncluster, $nearby_noncluster_points_count, $cluster_pts_near_noncluster_pt_count,
# 	  $max_d_to_consensus, $far_from_consensus_count, \%clusterpt_dist2consensus);
# }

# ##################

# sub least_noncluster_distance{	# call once for each cluster,
#   # to get least distance from cluster to outside cluster, etc.

#   my $thisclustids = shift;	# hash ref, keys ids in this cluster.
#   my $id1_id2ds = shift; # hash ref; key ids; value: array ref of strings with id2, distance12
#   my $near_noncluster_distance = shift; # set to some multiple of $link_max_distance
#   my $short_cluster_noncluster_edge_count = 0; # count number of cluster-noncluster edges < $near_noncluster_distance
#   my %near_noncluster_points = ();
#   my %cluster_pts_near_noncluster_pt = ();
#   my $min_distance_to_noncluster = 1;
#   #print STDERR "# in least... cluster size: ", scalar keys %$thisclustids, "\n";
#   for my $id1 (keys %$thisclustids) { # loop over elements of cluster
#     for my $s (@{$id1_id2ds->{$id1}}) { # loop over accessions with smallish distance to id1
#       my ($id2, $d) = @$s;
#       if ( !exists $thisclustids->{$id2} ) { # $id2 is not in this cluster
# 	if ($d < $min_distance_to_noncluster) {
# 	  $min_distance_to_noncluster = $d;
# 	}
# 	#	else{ # @$ is sorted by $d, small to large, so if this $d isn't small enough, we're done with @$s
# 	if ($d < $near_noncluster_distance) {
# 	  $short_cluster_noncluster_edge_count++;
# 	  $near_noncluster_points{$id2}++;
# 	  $cluster_pts_near_noncluster_pt{$id1}++;
# 	} else {
# 	  last;			#
# 	}
#       }
#     }
#   }
#   my $nearby_noncluster_points_count = scalar keys %near_noncluster_points;
#   my $cluster_pts_near_noncluster_pt_count = scalar keys %cluster_pts_near_noncluster_pt;
#   return ($min_distance_to_noncluster, $nearby_noncluster_points_count, $cluster_pts_near_noncluster_pt_count);
# }



# #####################
# sub cluster_quality1x{ # version which identifies 'dubious' ids/edges - so they can be pruned from clusters if desired.
#   my $accs = shift;    # array ref of ids of cluster members
#   my $edge_weight = shift; # hashref; keys: order id pairs, values: weights (distances)
#   my $graph = shift;
#   my $link_max_distance = shift;
#   my $output_line_string = '';
#   my @output_id_degree_maxds = ();
#   my @sorted_clusterids = sort {$a cmp $b} @$accs; # sort the accession ids in the cluster
#   my $cluster_size = scalar @sorted_clusterids;
#   my ($cluster_min_d, $cluster_max_d, $cluster_sum_d, $missing_distance_count, $intracluster_far_pair_count) = (10000, -1, 0, 0, 0);
#   my $sum_of_degrees = 0;
#   my $cluster_edge_count = 0;
#   my %id_maxd = map {$_ => -1} @sorted_clusterids;
#   my %id_sumdsq = map  {$_ => 0} @sorted_clusterids;
#   my $F = 0.2;
#   my $dubious_edges = {};
#   my %id_degree = ();
#   #  my %id_dists = map {$_ => []} @sorted_clusterids;
#   my $dubious_ids = {}; # keys: ids of cluster members with few links to other cluster members.
#   my $dubious;
#   while (my($i, $v) = each @sorted_clusterids) { # loop over every pair of ids in the cluster.
#     my $degree = $graph->degree($v); # degree of $v (=number of other nodes at dist <= $link_max_distance
#     $id_degree{$v} = $degree;
#     my $present_distance_count = 0;
#     if ($degree < $F*($cluster_size-1)) {
#       $dubious_ids->{$v} = 1;
#       $dubious = 1;
#     } else {
#       $dubious = 0;
#     }
#     my ($minw, $maxw) = (10000, -1);
    
#     for (my $j=$i+1; $j<scalar @sorted_clusterids; $j++) {
#       my $u = $sorted_clusterids[$j];
#       #  next if(exists $dubious_ids{$u});
#       if ($u eq $v) {
# 	next;
# 	warn "# $u $v  Why are they the same?\n";
# 	next;
#       }
#       my $edge_verts = ($v lt $u)? "$v $u" : "$u $v";
#       if ($dubious || $graph->degree($u) < $F*($cluster_size-1)) {
# 	printf("dubious: %s\n", $edge_verts) if(!exists $dubious_edges->{$edge_verts});
# 	$dubious_edges->{$edge_verts} = 1;
# 	$dubious_ids->{$u} = 1;
#       }
#       my $weight = $edge_weight->{$edge_verts} // -1;
#       if ($weight == -1) {
# 	$missing_distance_count++; # just counts pairs in cluster with distance not found
#       } else {
# 	$present_distance_count++;
# 	$minw = $weight if($weight < $minw);
# 	$maxw = $weight if($weight > $maxw);
# 	$id_maxd{$v} = max($id_maxd{$v}, $weight);
# 	$id_maxd{$u} = max($id_maxd{$u}, $weight);
# 	$id_sumdsq{$v} += $weight*$weight;
# 	$id_sumdsq{$u} += $weight*$weight;

# 	#	push @{$id_dists{$v}}, $weight;
# 	#	push @{$id_dists{$u}}, $weight;
# 	$cluster_sum_d += $weight;
# 	$cluster_edge_count++;
# 	$intracluster_far_pair_count++ if($weight > $link_max_distance); # count intra-cluster pairs with distance > $link_max_distance
#       }
#     }			 # loop over $u (only those with indices > $i)
#     $cluster_min_d = $minw if($minw < $cluster_min_d);
#     $cluster_max_d = $maxw if($maxw > $cluster_max_d);
#     $output_line_string .= "$v ";
#     # my $degree = $graph->degree($v);
#     $sum_of_degrees += $degree;
#     $output_line_string .= "$degree  "; # add number of edges joining cluster member to other cluster members.
#     #  printf("%d\n", scalar @{$id_dists{$v}});
#     #   push @output_id_degree_maxds, [$v, $degree, $id_maxd{$v}]; # $clusterid_dist2consensus->{$v}];
#     #}
#   } # loop over $v
#   while (my($id, $degree) = each %id_degree) {
    
#     my $rms_dist = sqrt($id_sumdsq{$id}/($cluster_size-1));
#     #    printf("id: %s  maxdist: %8.5f  rmsdist: %8.5f \n", $id, $id_maxd{$id}, $rms_dist);
#     push @output_id_degree_maxds, [$id, $degree, $id_maxd{$id}, $rms_dist];
#   }
  
#   @output_id_degree_maxds = sort {$b->[1] <=> $a->[1]} @output_id_degree_maxds;
#   return ($cluster_min_d, $cluster_sum_d/$cluster_edge_count, $cluster_max_d, $intracluster_far_pair_count, $missing_distance_count, \@output_id_degree_maxds, $dubious_edges, $dubious_ids);
# }

# sub auto_max_link_distance_xx{
#   #  my $edge_weight = shift;
#   my $id_closeidds = shift;
#   my $pow = shift;
#   my $f = shift;
#   my @distances = ();
#   while (my($id1, $id2ds) = each %$id_closeidds) {
#     push @distances, $id2ds->[0]->[1];
#     # for my $id2d (@$id2ds){
#     #   my ($id2, $d) = @$id2d;
#     #   push @distances, $d;
#     # }
#   }
#   print STDERR "## xx number of distances: ", scalar @distances, "\n";
#   my $cluster1d_obj = Cluster1d->new({label => '', xs => \@distances, pow => $pow});
#   print  "before two_cluster to choose cluster max distance\n";
#   print  "#  pow: $pow \n";
#   my ($Hopt, $maxQ, $Hopt_avg, $maxQ_avg, $Ledge, $Redge) = $cluster1d_obj->find_cluster_at_left($f);
#   my $Hmid_half_max = 0.5*($Ledge+$Redge);
#   my $H33 = (0.67*$Ledge + 0.33*$Redge);
#   # my ($Hoptx, $maxQx) = (0, 0); # $cluster1d_obj->two_cluster_x();
#   # $link_max_distance = $Hmid_half_max;
#   # print "# after two_cluster. Cluster max distance: $link_max_distance \n";
#   #  print "# before k-means/kde clustering \n";
#   #  my ($n_pts, $km_n_L, $km_n_R, $km_h_opt, $q, $kde_n_L, $kde_n_R, $kde_h_opt, $kde_q) = $cluster1d_obj->one_d_2cluster();
#   #  printf( STDERR "# clustering %5d points;  k-means: %5d below  %8.6f and  %5d above; q: %6.4f.  kde: %5d below  %8.6f  and %5d above; kde_q: %6.4f   Hopt: %6.4f  maxQ: %6.4f.  Hmhmx: %6.4f \n",
#   #       $n_pts, $km_n_L, $km_h_opt, $km_n_R, $q, $kde_n_L, $kde_h_opt, $kde_n_R, $kde_q, $Hopt, $maxQ, $Hmid_half_max);
#   printf(STDERR  "Hopt: %6.4f  maxQ: %6.4f. Hoptavg maxQavg: %6.4f %6.4f   Hmhmx: %6.4f H33: %6.4f \n", $Hopt, $maxQ, $Hopt_avg, $maxQ_avg, $Hmid_half_max, $H33);
#   return ($Hopt, $maxQ, $Ledge, $Redge, $Hopt_avg, $maxQ_avg);
# }
