#!/usr/bin/perl -w
use strict;

my $n_pedigrees = shift // 100;
my $self = shift // 0;

my $dsgs_in = 0;

my %accids = ();
while (<>) {		     # read accession ids from first 3 columns
  next if(/^\s*#/);
  my @cols = split(" ", $_);
  if ($cols[0] = 'MARKER') {
    $dsgs_in = 1;
    last;
  } elsif (scalar @cols == 3) {
    $accids{$cols[0]}++ if($cols[0] ne 'NA');
    $accids{$cols[1]}++ if($cols[1] ne 'NA');
    $accids{$cols[2]}++ if($cols[2] ne 'NA');
    last;
    print "X\n";

  }

print "Y\n";
  
}

while (<>) {
  next if(/^\s*#/);
  my @cols = split(" ", $_);
  if ($dsgs_in == 1) {
    $accids{$cols[0]}++ if($cols[0] ne 'NA');
  } elsif (scalar @cols == 3) {
    $accids{$cols[0]}++ if($cols[0] ne 'NA');
    $accids{$cols[1]}++ if($cols[1] ne 'NA');
    $accids{$cols[2]}++ if($cols[2] ne 'NA');
  } else {
    print STDERR "Unknown format. Columns: ", scalar @cols, "\n";
  }
}
my @acc_ids = keys %accids; 
my $n_accessions = scalar @acc_ids;
my $n_so_far =0;
while (1) {
  my $idx1 = int(rand($n_accessions));
  my $idx2 = int(rand($n_accessions));
  my $idx3 = ($self)? $idx2 : int(rand($n_accessions));
  next if($idx2 == $idx1  or  $idx3 == $idx1); # 'parents' are required be distinct from offspring.
  print  $acc_ids[$idx1], "   ", $acc_ids[$idx2], "  ", $acc_ids[$idx3], "\n";
  $n_so_far++;
  last if($n_so_far >= $n_pedigrees);
}
