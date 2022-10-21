#!/usr/bin/perl -w
use strict;

my $dosageA = shift // 0;
my $dosageB = shift // 1;
my $dosageC = shift // 3;

my $mdidx = 7;

my $n_accessions = 0;
my $first_line = <>;
my $line2 = <>;
my @line2_fields = split(" ", $line2);
my $line2_id = shift @line2_fields;
my $n_markers = scalar @line2_fields;

my @marker_dosage_counts; # = ([0, 0, 0, 0, 0, 0, 0, 0]) x $n_markers;

for (1..$n_markers) {
  push @marker_dosage_counts, [0, 0, 0, 0, 0, 0, 0, 0];
}
add_accession_gts_to_marker_dosage_counts(\@marker_dosage_counts, \@line2_fields);
$n_accessions++;

my @ids = ();
my %id_dosages = ();
push @ids, $line2_id;
$id_dosages{$line2_id} = \@line2_fields;

while (my $line = <>) {
  my @fields = split(" ", $line);
  my $id = shift @fields;
  add_accession_gts_to_marker_dosage_counts(\@marker_dosage_counts, \@fields);
  push @ids, $id;
  $id_dosages{$id} = \@fields;
  $n_accessions++;
}

# while (my($j, $v) = each @marker_dosage_counts) {
#   print STDERR "# ", join("  ", @$v), "\n";
# }
# exit;


# now check out all distinct triples of accessions:
while (my($i1, $par1id) = each @ids) {
  print STDERR "# $par1id \n";
  my $par1gts = $id_dosages{$par1id};
  for (my $i2 = $i1+1; $i2 < $n_accessions; $i2++) {
    my $par2id = $ids[$i2];
    # print STDERR "# $par1id $par2id \n";
    my $par2gts = $id_dosages{$par2id};
    my $markers_with_pattern = marker_indices_with_pattern($par1gts, $par2gts, $dosageA, $dosageB);
  #  print STDERR"#   $par1id  $par2id  ", scalar @$markers_with_pattern, "\n";
    for my $progid (@ids) {
      next if($progid eq $par1id  or  $progid eq $par2id);
      my $proggts = $id_dosages{$progid};
      my $count3 = 0;
      my $expect3 = 0;
      for my $index (@$markers_with_pattern, $dosageA, $dosageB) {
	my $dosage_counts = $marker_dosage_counts[$index];
	my $non_mds = $n_accessions - $dosage_counts->[$mdidx];
	if ($non_mds > 0) { # skip if all missing data for this marker
	  $expect3 += $dosage_counts->[$dosageC]/$non_mds;
	  $count3++ if($proggts->[$index] eq $dosageC);
	}
      }
      print "$par1id $par2id  $progid    $expect3  $count3  ", $count3/$expect3, "\n" if($expect3 > 0);
    }
  }
}




sub add_accession_gts_to_marker_dosage_counts{
  my $mrkrdoscounts = shift;
  my $dsgs = shift;
  while (my($i, $d) = each @$dsgs) {
    $d = $mdidx if($d eq 'NA');
    $mrkrdoscounts->[$i]->[$d]++;
  }
}

sub marker_indices_with_pattern{
  my $par1gts = shift;
  my $par2gts = shift;
  my $Adosage = shift;
  my $Bdosage = shift;
  my @pattern_indices = ();
  while (my($i, $v1) = each @$par1gts) {
    if ($v1 eq $Adosage) {
      my $v2 = $par2gts->[$i];
   #   print STDERR "#d1 d2:      $v1  $v2\n";
      if ($v2 eq $Bdosage) {
	push @pattern_indices, $i;
      }
    }
  }
  return \@pattern_indices;
}

