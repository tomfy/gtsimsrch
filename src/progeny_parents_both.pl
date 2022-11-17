#!/usr/bin/perl -w
use strict;
use Graph;

my $max_good_d = shift // 0.03;

# construct a graph. Nodes are accessions, edges point from parents to progeny.
# include only pedigrees with d <= $max_good_d

my %accids = ();		# keys accession ids, values: 1

my %acc_line = ();
my $g = Graph->new;
while (my $line = <>) {
  next if($line =~ /^\s*#/);
  my @cols = split(" ", $line);
  my ($id0, $id1, $id2, $d) = @cols[0,2,3,17];
  if ($d <= $max_good_d) {
    $acc_line{$id0} = $line;
    $accids{$id0} = 1; $accids{$id1} = 1; $accids{$id2} = 1;
    $g->add_edge($id1, $id0);	# edge from F parent to progeny
    $g->add_edge($id2, $id0);	# edge from M parent to progeny
  }
}

print scalar keys %acc_line, "\n";

for my $an_id (keys %acc_line) {
  my $n_in = $g->in_degree($an_id);
  my $n_out = $g->out_degree($an_id);
  if ($n_in eq 0) {
    if (exists $acc_line{$an_id}) {
      print $acc_line{$an_id};
    } else {
      #  print STDERR "$an_id  \n";
      # }
      print "$an_id  $n_in  $n_out \n";
    }
  }
}

while(my($an_id, $line) = each %acc_line){
  if(exists $accids{$an_id}){
    my $n_in = $g->in_degree($an_id);
    my $n_out = $g->out_degree($an_id);
    if($n_in eq 0){
      print $line;
    }
  }
}


