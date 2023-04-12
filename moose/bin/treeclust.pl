#!/usr/bin/perl -w
use strict;
use Graph::Undirected;
use Getopt::Long;
use List::Util qw(min max sum);
use File::Basename 'dirname';
use Cwd 'abs_path';
my ( $bindir, $libdir );
BEGIN {     # this has to go in Begin block so happens at compile time
  $bindir =
    dirname( abs_path(__FILE__) ) ; # the directory containing this script
  $libdir = $bindir . '/../lib';
  $libdir = abs_path($libdir);	# collapses the bin/../lib to just lib
}
print STDERR "[$bindir]  [$libdir] \n";
use lib $libdir;
use CNode;

my $max_d = shift // 0.15;

my %pair_distance = (); # keys: acc id pairs, e.g 'accession_a12 accession_z56'
my %forest = (); # keys: numerical tree ids; values: (root) node objects.
#my %accid_nodenumber = (); # keys: accession ids; values: numerical id of the leaf node representing the accession.
my %leaf = ();
my $node_number = 1;
my $factor = 2;

while(<>){
  next if(/^\s*#/);
  my @cols = split(" ", $_);
  my ($id1, $id2, $d) = @cols[0,1,5];
  #  print "$id1 $id2 $d\n";
  next if($d > $max_d);
  my $id_pair = ($id1 le $id2)? "$id1 $id2" : "$id2 $id1";
  my $new_node = CNode->new({id => $id1, node_number => $node_number,
			     Lchild => undef, Rchild => undef, parent => undef,
			    leaf_count => 1, distance => 0, p_distance => $max_d});
  $pair_distance{$id_pair} = $d;
  $forest{$node_number} = $new_node;
  $leaf{$id1} = $new_node;
  #$accid_nodenumber{$id1} = $node_number;
  $node_number++;
}

my @pairsdists = map([$_, $pair_distance{$_}], keys %pair_distance);
my @sorted_pairsdists = sort {$a->[1] <=> $b->[1]} @pairsdists;
# print $pairsdists[0]->[0], "  ", $pairsdists[0]->[1], "\n";
# print $sorted_pairsdists[0]->[0], "  ", $sorted_pairsdists[0]->[1], "\n";

# for my $pd (@sorted_pairsdists){
#   printf("%56s  %7.5f\n", @$pd);
# }

print "# number of d < $max_d pairs: ", scalar @sorted_pairsdists, "\n";
# sleep(1);

while(@sorted_pairsdists){
  my $next_link = shift @sorted_pairsdists;
  my ($id_pair, $d) = @$next_link;
  my ($id1, $id2) = split(" ", $id_pair);
  my $leaf1 = $leaf{$id1};
  my $leaf2 = $leaf{$id2};
  # my $Lnode_number = $accid_nodenumber{$id1};
  my $root1 = $leaf1->top_ancestor($leaf1); #
  #my $Rnode_number = $accid_nodenumber{$id2};
  my $root2 = $leaf2->top_ancestor($leaf2);
  # join the trees if the members of the closest pair are in distinct trees:
  if($root1->node_number() != $root2->node_number()){ # if the next closest pair belong to separate trees ...
    my $max12_d = max($root1->distance, $root2->distance);
    my $max12_size = max($root1->leaf_count(), $root2->leaf_count());
   print STDERR "$max12_d  $d  $max12_size  [" . $root1->id() . " ; " . $root2->id() . "]\n";
    if($max12_size == 1  or  $d < $factor*$max12_d){ # if at least one tree is not just a leaf, and
      print STDERR "            # join: ", $root1->id(), "  and  ", $root2->id(), "\n";
      my $new_node_id = $root1->id() . ' ' . $root2->id();
    #  my $new_parent_node = CNode->new({id => $new_node_id, node_number => $node_number});
    my $new_root = $root1->join($root2, $d, $node_number, $max_d); # create a parent node to root1 and root2, it is the root of the joined tree.
    $root1->p_distance($d);
    $root2->p_distance($d);
      #
    print STDERR "$node_number  ", $new_root->node_number(), "\n";
    $forest{$node_number} = $new_root; 
    delete $forest{$root1->node_number()}; # remove the 2 nodes which were the roots of
    delete $forest{$root2->node_number()}; # the trees which were just joined.
    $node_number++;
    }else{ # $d is big. Don't join, but update p_distance of $root1, $root2.
      $root1->p_distance($d);
      $root2->p_distance($d);
    }
  }
    
    # my $forest_size = scalar keys %forest;
    # sleep(2) if($node_number % 100  == 0);
    # print STDERR "forest size: $forest_size  node_number: $node_number \n";

}

while (my($n, $root) = each %forest) {
  if ($root->leaf_count() > 1) {
    if (1) {
      my $child_d = max($root->Lchild()->distance(), $root->Rchild()->distance());
      print $root->leaf_count(), "  ", $root->distance(), "  ",
	"  $child_d  " , $root->p_distance(), "  ", $root->id(), "\n";
    } else {			# output newick
      print $root->as_newick(), "\n";
    }
  }
}
