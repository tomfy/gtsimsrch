#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw(min max sum);

# vcf file should have one of :
# dosage DS,
# genotype GT, e.g. 0/0/0/1  or  0|0|0|1
# allele depth AD, e.g. 21,19
# and if it has GP (genotype prob.) or GQ (genotype quality) 
# we can also reject entries (i.e. regard as missing data) if these are not good enough

# usage:  vcf2dosage.pl <  <input vcf file>  >  <output file>

my $minGQ = 96;                 #
my $minGP = 0.9; # there must be 1 genotype with prob >= $minGP; i.e. one genotype must be strongly preferred.
my $transpose = 1; # default is to transpose; use -notrans to output untransposed.
# (duplicatesearch, find_parents require transposed)
# vcf: columns correspond to accessions, rows to markers

# if we don't believe can reliably resolve various heterozygous genotypes in polyploid case
# we can just lump together all heterozygous genotypes, map to just 3 genotypes:
my $map_to_012 = 0; # dosage = ploidy -> 2, 0 < dosage < ploidy -> 1, 0 -> 0, NA -> NA
my $field_to_use = 'AUTO'; # default is DS if present, then GT if present, then AD if present, then give up.
# recognized choices are DS (alt dosage e.g. 2), GT (genotype e.g. '/1/0' ) , AD (allele depths, e.g.'136:25' ), and AUTO.
my $ploidy = -1; # infer from data - user must specify if AD (allele depth) is specified.
my $hw = 0.33; # if not $map_to_012, round to integer if within +- $hw
my $delta = 0.1; # if map_to_012 [0, $delta ->0], [1-$hw, $ploidy-1+$hw] -> 1, [$ploidy-$delta, $ploidy] -> 2
my $min_read_depth = 1;

GetOptions(
	   'transpose!' => \$transpose, # -notranspose to output untransposed. (duplicatesearch requires transposed which is default)
	   'GQmin=f' => \$minGQ,	# min genotype quality.
	   'GPmin=f' => \$minGP,	# must
	   'field=s' => \$field_to_use,
	   'ploidy=f' => \$ploidy,
	   'map_to_012!' => \$map_to_012,
	   'hw|half_width=f' => \$hw,
	   'delta=f' => \$delta, 
	   'min_read_depth=f' => \$min_read_depth,
	  );

die "if specifying use of AD (allele depth), ploidy must also be specified\n" if($field_to_use eq 'AD'  and  $ploidy == -1);
#print STDERR "$transpose  $minGQ  $minGP  $field_to_use \n";

# read lines up to and including first starting with a single #
# that line has accession identifiers for columns 9, 10, ...
#
my @col_ids = ();
my @row_ids = ();
my @rows = ();
my @dosage_distribution = ();
while (<>) {
  my @cols = split(" ", $_);
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
# my $in_ploidy;
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
  push @row_ids, $row_id;	# store row (marker) id
  # print "$row_id  "; # the row id
  
  my @fields = split(':', $format_str);
  my $nfields = scalar @fields;
  my ($DSidx, $GTidx, $ADidx, $GPidx, $GQidx) = (-1, -1, -1, -1, -1);
  while (my($i, $f) = each @fields) {
    if ($f eq 'DS') {	# dosage. e.g. 0 or 1, but can be non-integer.
      $DSidx = $i;
    } elsif ($f eq 'GT') { # genotype. e.g. 0/1 (unphased) or 0|1 (phased), or 0/0/0/1 (tetraploid)
      $GTidx = $i;
    } elsif ($f eq 'AD') {	# allele depth. e.g. 142,31
      $ADidx = $i;
    } elsif ($f eq 'GP') { # genotype probability. e.g. 0.002,0.998,0.001   
      $GPidx = $i;
    } elsif ($f eq 'GQ') {	# genotype quality. e.g. 
      $GQidx = $i;
    }
  }
  @cols = @cols[9..$#cols]; # a typical elem: 0/0/0/0/0/1:79,18:97   GT:AD:DP  (genotype:allele depths:read depth)
  for my $e (@cols) {
    my $dosage = 'NA';		# indicates missing data
    my @field_values = split(":", $e);
    die if(scalar @field_values != $nfields);
    if ($field_to_use eq 'DS') {
      if ($DSidx >= 0) {      # DS; use dosage if present in vcf file.
	$dosage = $field_values[$DSidx];
      } else {
	die "Selected data field 'DS' not present.\n";
      }
    } elsif ($field_to_use eq 'GT') {
      if ($GTidx >= 0) {	# GT; use genotype if present
	if ($field_values[$GTidx] =~ /[.]/) { # no valid GT -> missing data
	} else {			      # count 0's and 1's
	  my $ref_allele_count = $field_values[$GTidx] =~ tr/0/0/;
	  my $alt_allele_count = $field_values[$GTidx] =~ tr/1/1/;
	  $dosage = $alt_allele_count;
	}
      } else {
	die "Selected data field 'GT' not present.\n";
      }
    } elsif ($field_to_use eq 'AD') {
      if ($ADidx >= 0) {
	# AD; use allele depth 
	my ($ref_depth, $alt_depth) = split(',', $field_values[$ADidx]);
	my $read_depth = $ref_depth + $alt_depth;
	if ($read_depth >= $min_read_depth) {
	  my $float_dosage = $ploidy*$alt_depth/$read_depth;
	  $dosage = int($float_dosage + 0.5);
	  if ($map_to_012) {
	    if (
		($float_dosage > $delta  and  $float_dosage < 1-$hw)
		or
		($float_dosage > $ploidy-1+$hw  and  $float_dosage < $ploidy-$delta)) {
	      $dosage = 'NA';
	    }
	  } else {
	    $dosage = 'NA' if(
			      abs($float_dosage - $dosage) > $hw
			      or
			      ($float_dosage > $delta  and  $float_dosage < 1-$hw)
			      or
			      ($float_dosage > $ploidy-1+$hw  and  $float_dosage < $ploidy-$delta)
			     );
	  }
	}
      } else {
	die "Selected data field 'AD' not present.\n";
      }
    } else {   # no data field selected, look for DS, then GT, then AD
      if ($DSidx >= 0) {      # DS; use dosage if present in vcf file.
	$dosage = $field_values[$DSidx];
      } elsif ($GTidx >= 0) {	# GT; use genotype if present
	if ($field_values[$GTidx] =~ /[.]/) { # no valid GT -> missing data
	} else {			      # count 0's and 1's
	  my $ref_allele_count = $field_values[$GTidx] =~ tr/0/0/;
	  my $alt_allele_count = $field_values[$GTidx] =~ tr/1/1/;
	  $dosage = $alt_allele_count;
	}
      } elsif ($ADidx >= 0) {	# AD; use allele depth 
	my ($ref_depth, $alt_depth) = split(',', $field_values[$ADidx]);
	my $read_depth = $ref_depth + $alt_depth;
	#	$dosage = int($ploidy*$alt_depth/$read_depth + 0.5) if($read_depth > 0);
	if ($read_depth >= $min_read_depth) {
	  my $float_dosage = $ploidy*$alt_depth/$read_depth;
	  $dosage = int($float_dosage + 0.5);
	  $dosage = 'NA' if(abs($float_dosage - $dosage) > $hw);
	}
      }
    }
  

    ###   do quality checks here if GP or GQ available   ###
    # check if GQ present but too low:
    $dosage = 'NA' if($GQidx >= 0  and  $field_values[$GQidx] < $minGQ);
    # check if GP present but insufficiently strong preference for one genotype:
    if ($GPidx >= 0) {
      my @gt_probs = split(',', $field_values[$GPidx]);
      $dosage = 'NA' if(max(@gt_probs) < $minGP);
    }
    if ($dosage eq 'NA') {
      $missing_data_count++ 
    } else {			# $dosage ne 'NA'
      $ploidy = $dosage if($dosage > $ploidy);
      $dosage_distribution[$dosage]++;
    }
   

    push @dosages_this_row, $dosage;
  }				# end loop over entries in a row
  die "# n col ids: ", scalar @col_ids, "; n dosages in row: ", scalar @dosages_this_row, "\n" if(scalar @dosages_this_row != scalar @col_ids);
  #print "# n col ids: ", scalar @col_ids, "; n dosages in row: ", scalar @dosages_this_row, "\n";
  #sleep(1);
  push @rows, \@dosages_this_row;

}				# end loop over rows

# #####  output  #####
my ($transtr, $n_rows_out, $n_cols_out) = ($transpose)?
  ("# transpose", scalar @col_ids, scalar @row_ids) : ("# no transpose", scalar @row_ids, scalar @col_ids);
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
      my $d = $r->[$i];
      # at this point $d = 0, 1, 2, ... , ploidy or NA (missing data)
      if ($map_to_012) {
	if ($d eq 'NA') {
	  # no change
	} elsif ($d eq $ploidy) {
	  $d = 2;
	} elsif ($d > 0) {
	  $d = 1;
	}			# else $d = 0, no change
      }
      print "$d ";
    }
    print "\n";
  }
}
