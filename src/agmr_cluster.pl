#!/usr/bin/perl -w
use strict;
use Graph::Undirected;


# read simsearch output and find clusters of accessions
# with near-identical genotypes
# specifically, create a graph, and for any pairs in simsearch output with
# agmr <= $max_agmr, make an edge in the graph between the accessions of the pair.
# the clusters are then the connected components of the graph.
# i.e. two accession belong to the same cluster iff you can get from
# one to the other by traversing edges of the graph.
# one line of output for each cluster, with fields
# cluster size, min agmr, max agmr, 'nbad (should be 0)', the ids of accessions in cluster.
# e.g.:    3  0.0104 0.0142 0  TMEB778:250254008  TMEB778:250304613  TMEB778:250597946
# the simsearch output should have all the agmrs for all the pairs in each cluster,
# and in that case the 4th field ('nbad') will be 0. If nbad > 0,
# consider rerunning this script with a smaller $max_agmr,
# or rerun simsearch with a larger value of 'max estimated agmr'

# usage:
# agmr_cluster.pl < simsrch.out  >  agmr_clusters

my $max_agmr = shift // 0.04; # default clustering parameter - max edge weight of graph
my $g = Graph::Undirected->new;

my %edge_weight = ();		# keys: vertex pairs, values: agmrs
while (my $line = <>) {
  next if($line =~ /^\s*#/);
  my ($id1, $id2, $usable_chunks, $match_chunks, $est_agmr, $agmr) = split(" ", $line);
  my $edge_verts = ($id1 lt $id2)? "$id1 $id2" : "$id2 $id1";
  $edge_weight{$edge_verts} = $agmr;
  if ($agmr <= $max_agmr) {
    $g->add_weighted_edge($id1, $id2, $agmr);
  }
}

my @ccs = $g->connected_components;

my @output_lines = ();

my $count = 0;
for my $acc (@ccs) {
  my $cc_size = scalar @$acc;
  $count += $cc_size;
  my $output_line_string = '';
  my @sorted_cc = sort {$a cmp $b} @$acc; # sort the accession ids in the cluster
  my ($ccminw, $ccmaxw, $nbad) = (100000, -1, 0);
    while(my($i, $v) = each @sorted_cc){
    my ($minw, $maxw) = (100000, -1);
      for(my $j=$i+1; $j<scalar @sorted_cc; $j++){
	my $u = $sorted_cc[$j];
      next if($u eq $v);
      my $edge_verts = ($v lt $u)? "$v $u" : "$u $v";
      my $weight = $edge_weight{$edge_verts} // -1;
      if ($weight == -1) {
	$nbad++;
      } else {
	$minw = $weight if($weight < $minw);
	$maxw = $weight if($weight > $maxw);
      }
    }
    $ccminw = $minw if($minw < $ccminw);
    $ccmaxw = $maxw if($maxw > $ccmaxw);
    $output_line_string .= "$v  ";
  }
  $output_line_string = "$cc_size  $ccminw $ccmaxw $nbad  " . $output_line_string . "\n";
  push @output_lines, $output_line_string;
}

my @sorted_output_lines = sort { compare_str($a, $b) }  @output_lines;
print "# graph max edge length: $max_agmr. Found ", scalar @ccs, " groups, with total of $count accessions.\n";
print join('', @sorted_output_lines);

sub compare_str{ # sort by size of cluster, tiebreaker is id of first accession in cluster. 
  my $str1 = shift;
  my $str2 = shift;
  my @cols1 = split(" ", $str1);
  my @cols2 = split(" ", $str2);
  return ($cols1[0] != $cols2[0])? $cols1[0] <=> $cols2[0] : $cols1[4] cmp $cols2[4];
}
