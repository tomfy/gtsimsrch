#!/usr/bin/perl -w
use strict;
use Graph::Directed;

# construct DAG from pedigree_test output
# output ids of triples with selected relationship:
# 1: accession, parent1, parent2,
# 2: accession, progeny1, progeny2
# 3: accession, parent, progeny

# typical usage:  triple_3_cases.pl  2 < simsrch.out

my $case = shift // 1;

# initialize graph:

my $g = Graph::Directed->new;
my %ids = (); # just store all the ids of the nodes in the graph

while(<>){ # read pedigree_test output file and add edges to the graph
    next if(/^\s*#/); # skip comment lines
    my @cols = split(" ", $_);
    my ($aid, $pid1, $pid2) = @cols[0,2,3];
    $g->add_edge($pid1, $aid);
    $g->add_edge($pid2, $aid);
   
    $ids{$aid} = 1;
    $ids{$pid1} = 1;
    $ids{$pid2} = 1;
}

# check if it is a DAG:
#my $has_a_cycle = $g->has_a_cycle();
#my @cycle = $g->find_a_cycle(); print join(" ", @cycle), "\n";

my $is_DAG = $g->is_dag();
print "# the graph is ", ($is_DAG)? "a DAG" : "NOT a DAG", "\n";

for my $anid (keys %ids){
    my $in_degree = $g->in_degree($anid);
    die "$anid has in_degree of $in_degree\n" if($in_degree > 2);
    my $out_degree = $g->out_degree($anid);
    if($case == 1){ 
	if($in_degree == 2){
	    my $str = "$anid ";
	    my @edges_to = $g->edges_to($anid); # list of edges to the node;
	    # edge is array ref of 2 node ids; directed edge goes from zeroth to first element
	    for my $e (@edges_to){
		$str .= " " . $e->[0];
		# $str .= join(" ", @$e) . " ";
	    }
	    print "$str \n";
	}else{
	    # too few parents known, can't do case 1
	}
    }elsif($case == 2){
	if($out_degree >= 2){
	    # my $str = "$anid ";
	    my @edges_from = $g->edges_from($anid);
	    #for my $e (@edges_from){
	    for(my $i=0; $i < $#edges_from; $i+=2){
		my $e1 = $edges_from[$i];
		my $e2 = $edges_from[$i+1];
		my $str = "$anid " . $e1->[1] . " " . $e2->[1];
		print "$str\n";
	    }
	}
    }elsif($case == 3){
	if($in_degree >=1  and $out_degree  >= 2){
	    my @edges_to = $g->edges_to($anid);
	    my @edges_from = $g->edges_from($anid);
	    my $in_edge = $edges_to[0];
	    for my $out_edge (@edges_from){
		my $str =  "$anid " . $in_edge->[0] . " " . $out_edge->[1];
		print "$str\n";
	    }
	}
    }else{
	die "unknown case: $case\n"
    }

}



    
    