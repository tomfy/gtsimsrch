#!/usr/bin/perl -w
use strict;

# input a vcf file and a list of accession ids to select,
# outputs a vcf file with just those accessions (columns).

# usage:  perl vcf_acc_subset.pl   <ids_to_keep_file>  <  <vcf file>  >  <output file>

my $ids_file = shift;

my %ids = ();
open my $fhids, "<", "$ids_file";
while (<$fhids>) {
  if (/(\S+)/) {
    $ids{$1}++;
  }
}
my @ids_to_keep = keys %ids;
# print scalar @ids_to_keep, "\n";

my @accids;
my %id_col = ();
# my %col_id = ();
while (<>) { # read in vcf file
  #print;
  next if(/^\s*##/);
  if (/^\s*#/) {
    print STDERR "# acc id line.\n";
    @accids = split("\t", $_);
    my $init_stuff = join("\t", @accids[0..8]);
    @accids = @accids[9..$#accids]; # keep just the accids
    print "$init_stuff";
    last;
  }
}
my $first_accid = $accids[0];
print STDERR scalar @accids, "  ", $accids[0], "\n";

while (my($i, $accid) = each @accids) {
  $id_col{$accid} = $i;		# store the col of each accession id.
 #  print STDERR "### $accid  $i\n";
  # $col_id{$i} = $accid;
}
#sleep(10);
for my $anid (@ids_to_keep) {
  
  my $acol = $id_col{$anid};
#  print STDERR "$anid  \n";
#  print STDERR "  $acol\n";
#  sleep(1);
  print "\t", $accids[$acol];
}
print "\n"; # exit;

print STDERR "first accid: $first_accid    col: ", $id_col{$first_accid}, "\n";

print STDERR "# ", scalar keys %id_col, "  \n"; # , scalar keys %col_id, "\n";
#exit;

my $marker_line_count = 0;
while (<>) {
  my @gts = split("\t", $_); # the genotype information for all acessions, one marker.
  my $init_stuff = join("\t", @gts[0..8]);
  print "$init_stuff";
  @gts = @gts[9..$#gts];
  for my $anid (@ids_to_keep) {
    print STDERR "id to keep $anid\n";
   # sleep(1);
    my $acol = $id_col{$anid};
    print "\t", $gts[$acol];
  }
  print "\n";
}
