#!/usr/bin/perl -w
use strict;

# read in a dosage matrix file
# and report some per accession info,
# as well as summary info on whole set
# number of accessions, markers,
# frequencies of each dosage and missing data

my $n_marker_ids;
my $n_accessions = 0;
my ($cume_zero_count, $cume_one_count, $cume_two_count, $cume_X_count, $cumulative_bad_accession_count) = (0, 0, 0, 0, 0);
print "# dosage counts for each accession.\n";   
print "# acc_no                             acc_id        n0      n1      n2      nX   n_markers     maf    heterozygosity\n";
while (<>) {
  next if(/^\s*#/);
  if (/^MARKER/) {
    my @cols = split(" ", $_);
    # while(my($ii, $xx) = each @cols){
    #   print STDERR "###   $xx  $ii  $ii \n";
    # }
    $n_marker_ids = scalar @cols - 1;
    print "# n marker ids: $n_marker_ids\n";
  }elsif(/^CHROMOSOME/) { # ignore this line
    
  } else {
    s/^\s*(\S+)//; # delete first field (accession id)
   my $acc_id = $1;
 #  while (/\s+\S/g) { $marker_count++ } # slow compared to using tr as in following 4 lines.
    my $zero_count = $_ =~ tr/0/0/;
    my $one_count = $_ =~ tr/1/1/;
    my $minus_count = $_ =~ tr/-/-/;
    my $plus_count = $_ =~ tr/+/+/;
    my $two_count = $_ =~ tr/2/2/;
    my $X_count = $_ =~ tr/X/X/;
    #my $acc_heterozygosity = $one_count/($zero_count + $one_count + $two_count);
    # print STDERR "$acc_id   $zero_count  $one_count ($minus_count, $plus_count)  $two_count  $X_count  $acc_heterozygosity\n";
    my $marker_count =  $zero_count+$one_count+$two_count+$X_count;
    die "# Marker count discrepancy: n_marker_ids $n_marker_ids  n_markers: $marker_count\n" if($marker_count != $n_marker_ids);
    $cume_zero_count += $zero_count;
    $cume_one_count += $one_count;
    $cume_two_count += $two_count;
    $cume_X_count += $X_count;
    $n_accessions++;
    my $ok_count = $marker_count - $X_count;
    $cumulative_bad_accession_count++ if($ok_count == 0);
    my $alt_allele_freq = ($ok_count > 0)? ($one_count + 2*$two_count)/(2*$ok_count) : -1;
    my $maf = ($alt_allele_freq <= 0.5)? $alt_allele_freq : 1.0 - $alt_allele_freq;
    my $heterozygosity = ($ok_count > 0)? $one_count/$ok_count : -1;
    printf(STDOUT "%5ld %30s   %6ld %6ld (%5d,%5d) %5ld %7ld   %7ld   %8.4f    %8.4f\n",
	   $n_accessions, $acc_id, $zero_count, $one_count, $minus_count, $plus_count, $two_count, $X_count, $marker_count, $maf, $heterozygosity);
   # print STDOUT "#           ", $one_count + 2*$two_count, "  ", $ok_count, "  ", $alt_allele_freq, "\n"
  }
}
print "# n_accessions: $n_accessions  n_markers: $n_marker_ids\n";
my $n_total_gts = $n_accessions*$n_marker_ids;
print "# dosage frequencies:\n";
print "#    zero     one     two       X\n";
printf(STDOUT "#  %6.4f  %6.4f  %6.4f  %6.4f\n",
  $cume_zero_count/$n_total_gts, $cume_one_count/$n_total_gts,
       $cume_two_count/$n_total_gts, $cume_X_count/$n_total_gts);
print STDOUT "# accessions with all missing data: $cumulative_bad_accession_count\n";
