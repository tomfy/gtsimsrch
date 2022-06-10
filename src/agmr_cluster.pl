#!/usr/bin/perl -w
use strict;
use Graph::Undirected;

my $max_agmr = shift // 0.05;
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
  my @sorted_cc = sort @$acc;
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

my @sorted_output_lines = sort {get_first_number($a) <=> get_first_number($b)} @output_lines;
print "# found ", scalar @ccs, " groups, with total of $count accessions.\n";
print join('', @sorted_output_lines);

sub get_first_number{
  my $str = shift;
  $str =~ /^(\d+)/;
  my $result = $1 // -1;
  return $result;
}
