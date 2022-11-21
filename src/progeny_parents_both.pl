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
  if ($d <= $max_good_d) {	# pedigree looks good, add to graph
    $acc_line{$id0} = $line;
    $accids{$id0} = 1; $accids{$id1} = 1; $accids{$id2} = 1;
    $g->add_edge($id1, $id0);	# edge from F parent to progeny
    $g->add_edge($id2, $id0);	# edge from M parent to progeny
  }
}

#qprint scalar keys %acc_line, "\n";

# for my $an_id (keys %accids) {
#   my $n_in = $g->in_degree($an_id);
#   my $n_out = $g->out_degree($an_id);
#  # print STDERR "n in: $n_in \n";
#   if ($n_out eq 0) {
#     if (exists $acc_line{$an_id}) {
#       print $acc_line{$an_id} // "XXXX\n";
#     } else {
#       #  print STDERR "$an_id  \n";
#       # }
#       print "$an_id  $n_in  $n_out \n";
#     }
#   }
# }

open my $fhprogonly, ">", "progeny_only";
open my $fhboth, ">", "both";
open my $fhparentsonly, ">", "parents_only";

for my $an_id (keys %accids) {
  my $n_in = $g->in_degree($an_id);
  my $n_out = $g->out_degree($an_id);
  if (exists $acc_line{$an_id}) {
    my $line = $acc_line{$an_id};
    chomp $line;
    if ($n_out > 0) {		# has progeny
      if ($n_in > 0) {		# has parent
	print $fhboth $line . "  $n_in $n_out\n";
      } else {
	print $fhprogonly $line . "  $n_in $n_out\n";
      }
    } else {			# no progeny
      if ($n_in > 0) {
	print $fhparentsonly $line . "  $n_in $n_out\n";
      } else {
	die "acc: $an_id has neither parents nor progeny?\n";
      }
    }
  } else {
    print $fhprogonly "$an_id $n_in $n_out\n";
  }
}

