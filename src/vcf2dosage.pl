#!/usr/bin/perl -w
use strict;
use Getopt::Long;

# vcf file should have one of :
# dosage DS, 
# genotype GT, e.g. 0/0/0/1  or  0|0|0|1
# allele depth AD, e.g. 21,19
# and if it has GP (genotype prob.) or GQ (genotype quality) 
# we can also reject entries if these are not good enough (i.e. regard as missing data)

#my $out_ploidy = shift // -1; # default: same as in_ploidy, which will be auto-detected

my $minGQ = 96;                 #
my $maxGP = 0.1; # GP required to be <= $maxGP or >= 1-$maxGP , i.e. one allele strongly preferred.
my $transpose = 1; # default is to transpose; use -notrans to output untransposed.

my $ploidy = -1;

GetOptions(
	   'transpose!' => \$transpose, # -noplot to suppress plot - just see histogram as text.
	   'GQmin|qmin' => \$minGQ,	# min genotype quality.
	   'GPmax|delta' => \$maxGP, # either ref. or alt. allele required to have prob. less than this.
	  );



# read lines up to and including first starting with a single #
# that lines has identifiers for columns 9, 10, ...
#
my @col_ids = ();
my @row_ids = ();
my @rows = ();
my @dosage_distribution = ();
while (<>) {
  next if(/^\s*##/);
  if (/^\s*#/) {
    @col_ids = split(" ", $_);
    @col_ids = @col_ids[9..$#col_ids];
    #  print "MARKER  ", join(' ', @col_ids);
    last;
  }
}

#print "# number of col ids: #  ", scalar @col_ids, "\n";
#sleep(2);
# read the rows with genotype data
my $in_ploidy;
my $format_string = 'xxx';
my $missing_data_count = 0;
while (<>) {
  my @dosages_this_row = ();
  my @cols = split(" ", $_);
  my $row_id = $cols[2];
  my $format_str = $cols[8];
  # if($format_str != $format_string){
  #   if($format_string eq 'xxx'){
  #     $format_string = $format_str;
  #   }else{
  #    die;
  #   }
  # }
  push @row_ids, $row_id;
  # print "$row_id  "; # the row id
  
  my @fields = split(':', $format_str);
  my $nfields = scalar @fields;
  my ($DSidx, $GTidx, $ADidx, $GPidx, $GQidx) = (-1, -1, -1, -1, -1);
  while (my($i, $f) = each @fields) {
    if ($f eq 'DS') {
      $DSidx = $i;
    } elsif ($f eq 'GT') {
      $GTidx = $i;
    } elsif ($f eq 'AD') {
      $ADidx = $i;
    } elsif ($f eq 'GP') {
      $GPidx = $i;
    } elsif ($f eq 'GQ') {
      $GQidx = $i;
    }
  }
  @cols = @cols[9..$#cols]; # a typical elem: 0/0/0/0/0/1:79,18:97   GT:AD:DP  (genotype:allele depths:read depth)
  for my $e (@cols) { 
    my $dosage = 'NA';
    my @field_values = split(":", $e);
    die if(scalar @field_values != $nfields);
    if ($DSidx >= 0) {		# use DS if present in vcf file.
      $dosage = $field_values[$DSidx];
    } elsif ($GTidx >= 0) {			# use GT if present
      if (! ($field_values[$GTidx] =~ /[.]/)) { # if not '.' present if GT
	my $refallele_count = $field_values[$GTidx] =~ tr/0/0/;
	my $altallele_count = $field_values[$GTidx] =~ tr/1/1/;
	$dosage = $altallele_count;
      }
    } elsif ($ADidx >= 0) {
      my ($ref_depth, $alt_depth) = split(',', $field_values[$ADidx]);
      my $read_depth = $ref_depth + $alt_depth;
      $dosage = int($alt_depth/$read_depth + 0.5) if($read_depth > 0);
    }
    # check if GQ present but too low:
    $dosage = 'NA' if($GQidx >= 0  and  $field_values[$GQidx] < $minGQ);
    # check if GP present but insufficiently strong preference for one allele:
    if ($GPidx >= 0) {
      my ($ref_prob, $alt_prob) = split(',', $field_values[$GPidx]);
      $dosage = 'NA' if($ref_prob > $maxGP  and  $alt_prob > $maxGP);
    }
    if ($dosage eq 'NA') {
      $missing_data_count++ 
    } else {			# $dosage ne 'NA'
      $ploidy = $dosage if($dosage > $ploidy);
      $dosage_distribution[$dosage]++;
    }
    ###   do quality checks here if GP or GQ available   ###

    push @dosages_this_row, $dosage;
  }				# end loop over entries in a row
  die "# n col ids: ", scalar @col_ids, "; n dosages in row: ", scalar @dosages_this_row, "\n" if(scalar @dosages_this_row != scalar @col_ids);
  #print "# n col ids: ", scalar @col_ids, "; n dosages in row: ", scalar @dosages_this_row, "\n";
  #sleep(1);
  push @rows, \@dosages_this_row;

}				# end loop over rows

# #####  output  #####
my ($transtr, $n_rows_out, $n_cols_out) = ($transpose)? ("# transpose", scalar @row_ids, scalar @col_ids) : ("# no transpose", scalar @col_ids, scalar @row_ids);
my $n_elements = $n_rows_out * $n_cols_out;
my $info_string = "$transtr\n";
$info_string .= sprintf("# outputting %d rows and %d columns of data.\n", $n_rows_out, $n_cols_out);
$info_string .= sprintf("# ploidy: %d\n", $ploidy);
$info_string .= sprintf("# dosage distribution:\n");
for my $i (0..$ploidy) {
  $info_string .= sprintf("#    %2d        %8d   %8.6f\n", $i, $dosage_distribution[$i], $dosage_distribution[$i]/$n_elements);
}
$info_string .= sprintf("# missing data %8d   %8.6f\n", $missing_data_count, $missing_data_count/$n_elements);
print STDERR $info_string;
print STDOUT $info_string;
if (! $transpose) {
  print "MARKER ", join(" ", @col_ids), "\n";
  die if(scalar @row_ids  !=  scalar @rows);
  while (my ($i, $row_dosages) = each(@rows)) {
    print $row_ids[$i], "  ", join(" ", @$row_dosages), "\n";
  }
} else {			# transpose
  print "MARKER ", join(" ", @row_ids), "\n";
  while (my ($i, $col_id) = each @col_ids) {
    # print col i as a row
    print "$col_id  ";
    while (my($j, $r) = each @rows) {
      print $r->[$i], " ";
    }
    print "\n";
  }
}
