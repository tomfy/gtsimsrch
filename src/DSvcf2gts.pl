#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw(min max sum);
use File::Spec 'splitpath';
use Devel::Size qw(size total_size);

# This version is for a vcf file with dosages only
# these dosages are accepted if within $delta of an integer
# otherwise it is considered to be missing data.

#

# genotype GT, e.g. 0/1  or  0|1
# dosage DS,
# allele depth AD, e.g. 21,19
# and if it has GP (genotype prob., e.g. 0.9,0.16,0.04) or GQ (genotype quality, e.g. 98) 
# we can also reject entries (i.e. regard as missing data) if these are not good enough

# usage:  vcf2gts.pl -i <input vcf file>  -o <output file>

# vcf: columns correspond to accessions, rows to markers

# if we don't believe can reliably resolve various heterozygous genotypes in polyploid case
# we can just lump together all heterozygous genotypes, map to just 3 genotypes:
#my $map_to_012 = 0; # dosage = ploidy -> 2, 0 < dosage < ploidy -> 1, 0 -> 0, X -> X
my $field_to_use = 'DS'; # default GT. If requested field is not present exit (or ???)
# recognized choices are  GT (genotype e.g. '/1/0' ), DS (alternative allele dosage e.g. 2), AD (allele depths, e.g.'136:25' ).
my $ploidy = 2; # need to specify, will die if find dosages greater than specified ploidy.-
my $inferred_ploidy = -1; # infer from data; die if > specified $ploidy
my $delta = 0.1; # if not $map_to_012, round to integer if within +- $delta
                 # if map_to_012 [0, $delta ->0], [1-$delta, $ploidy-1+$delta] -> 1, [$ploidy-$delta, $ploidy] -> 2
my $input_vcf_filename = undef;
my $output_genotypes_filename = undef; # default: construct from input filename
my $missing_data_string = 'X';
#my $minGQ = 0;	 # if GQ present, must be >= this.
#my $minGP = 0.0; # if GP present, there must be 1 genotype with prob >= $minGP; i.e. one genotype must be strongly preferred.
my $min_marker_avg_pref_gt_prob = -1.0; # default is negative (meaning don't filter on this)
my $max_marker_missing_data_fraction = 1.0; # remove markers with excessive missing data. Default is keep all.
my $min_marker_maf = 0;
my $info_string = "# command: " . join(" ", @ARGV) . "\n";

GetOptions(
	   'input_file|vcf=s' => \$input_vcf_filename,
	   'output_file|dosage_file=s' => \$output_genotypes_filename,
	  # 'GQmin=f' => \$minGQ,	# min genotype quality.
	  # 'GPmin=f' => \$minGP,	# 
	   'delta=f' => \$delta, 
	  # 'min_read_depth=f' => \$min_read_depth,
	   'max_marker_md_fraction=f' => \$max_marker_missing_data_fraction,
	   'min_maf|min_marker_maf=f' => \$min_marker_maf,
	   'ploidy=f' => \$ploidy,
	  # 'map_to_012!' => \$map_to_012,
	  );

# #####  check for input filename; construct output filename if not specified  #####
# die "if specifying use of AD (allele depth), ploidy must also be specified\n" if($field_to_use eq 'AD'  and  $ploidy == -1);
die "Must specify input vcf filename. \n" if(!defined $input_vcf_filename);

if (!defined $output_genotypes_filename) { # construct an output filename from input vcf file name.
  (my $vol, my $dir, $output_genotypes_filename) = File::Spec->splitpath($input_vcf_filename);
  if ($output_genotypes_filename =~ /vcf$/) {
    $output_genotypes_filename =~ s/vcf/dosages/; #  replace final 'vcf' with 'dosages'
  } else {			# just append '.dosages'
    $output_genotypes_filename .= ".dosages";
  }
}

my $output_markerids_filename;

open my $fhout, ">", "$output_genotypes_filename" or die "Couldn't open $output_genotypes_filename for writing.\n";
# ##########################################################

# #####  output run parameters  #####

$info_string .= "# output file: $output_genotypes_filename\n";
$info_string .= "# data field to use: $field_to_use\n";
$info_string .= "# delta: $delta\n";
#$info_string .= "# GPmin: $minGP ; GQmin: $minGQ\n";
$info_string .= "# max marker missing data fraction: $max_marker_missing_data_fraction \n";
# $info_string .= "\n";
print $fhout $info_string;
print STDERR $info_string;
$info_string = '';
# ##################################################

# #####  read in initial comments and marker ids  #####
# read lines up to and including first starting with a single  '#'
# that line has accession identifiers in columns 9, 10, ...

my @col_ids = ();
my @dosage_distribution = ();
open my $fhin, "<", "$input_vcf_filename" or die "Couldn't open $input_vcf_filename for reading.\n";
while (<$fhin>) {
  next if(/^\s*##/);
  s/\s+$//;
  if (/^\s*#/) {
    @col_ids = split("\t", $_);
    my $n_col_ids = scalar @col_ids;

    @col_ids = @col_ids[9..$#col_ids];
    while (my($i, $v) = each @col_ids) {
      $col_ids[$i] //= '-';
      $col_ids[$i] =~ s/\s+/_/; # replace whitespace with underscore
    }
    last;
  }
}
print $fhout "# number of accession ids: #  ", scalar @col_ids, "\n"; # if(! $plink_format);
print STDERR "# number of accession ids: #  ", scalar @col_ids, "\n";

# #####  done reading accession ids  #####


# #####  read the rows with genotype data  #####
my @row_ids = ();		# i.e. marker ids
my @alt_row_ids = ();
my %rowid = ();
my %altrowid = ();
# my @rows = ();
#my @row_strings = ();
my $missing_data_count = 0;
my $markers_read_count = 0;
my @accession_gt_strings = (); # ith element with be string of genotypes for the ith accession
my %fieldused_count = (); # keys: 'DS', 'GT', etc.; values: count of gt data of that format.
my $n_accessions = 0;
while (<$fhin>) {
  my $marker_missing_data_count = 0;
  my $marker_alt_allele_count = 0;
  my @marker_dosage_distribution = ();
  my @genotypes_this_row = ();

  s/\s*$//; # remove any whitespace at end.
  my @cols = split("\t", $_);

  my $row_id = $cols[2]; # store row id in file (if these are not distinct, we will use $alt_row_id).
  my $alt_row_id = $cols[0] . "_" . $cols[1]; # construct a marker id from col[0] (chromosome number) and col[1] (position)
  my $chr_number = $cols[0];
  $chr_number =~ s/^chr//;

  my $format_str = $cols[8]; # e.g. 'GT:DS:AD' tells which types of data are present.
  my ($allele_zero, $allele_one) = @cols[3,4]; # e.g. 'A' and 'T' 

  my ($avg_pref_gt_prob, $pgtprob_count) = (0, 0); 

  # record which data fields are present (DS, GT, AD, GP, GQ)
  die "format string must be 'DS' is $format_str. Exiting\n" if($format_str ne 'DS'); 
 # my @fields = split(':', $format_str);
 # print "fields: ", join(", ", @fields), "\n";
 # my $nfields = scalar @fields;
 # my ($DSidx, $GTidx, $ADidx, $GPidx, $GQidx) = fields_present($format_str);

  # for this row (i.e. marker) loop over data for all accessions. Typical elements:
  # (GT:AD:DP): 0/0/0/0/0/1:79,18:97   ((hexaploid)genotype:allele depths:read depth)
  # (GT:DS:GP): 0|0:0.001:0.999,0.001,0    ((diploid)genotype:dosage:est_genotype_prob)
  @cols = @cols[9..$#cols];	# cols 9 and above have genotype info.
  # for my $e (@cols) {
  $n_accessions = 0;
  my @dosages_this_marker = ();
  while (my($acc_idx, $e) = each @cols) {
    my $dosage = $missing_data_string; # indicates missing data
	if (length $e > 0) {
	  my $float_dosage = $e;
	  my $int_dosage = int($float_dosage + 0.5); # round to integer
	  if (abs($float_dosage - $int_dosage) <= $delta) { # float_dosage is close to integer, use
	    $dosage = $int_dosage;
	  }			# else regard as missing data.
	}

    if ($dosage eq $missing_data_string) {
      $marker_missing_data_count++;
    } else {			# $dosage ne $missing_data_string
      $ploidy = $dosage if($dosage > $ploidy);
      $marker_dosage_distribution[$dosage]++;
      $marker_alt_allele_count += $dosage;
    }
    push @dosages_this_marker, $dosage;
    $n_accessions++
  } 	# end loop over entries in a row

  die "Inferred ploidy ($inferred_ploidy) is greater than specified ploidy ($ploidy).\n" if($inferred_ploidy > $ploidy);
  my $marker_total_allele_count = $ploidy*($n_accessions - $marker_missing_data_count);
  my $marker_alt_allele_frequency = $marker_alt_allele_count/$marker_total_allele_count;
  my $marker_minor_allele_frequency = ($marker_alt_allele_frequency <= 0.5)?
    $marker_alt_allele_frequency :
    1.0 - $marker_alt_allele_frequency;

  $missing_data_count += $marker_missing_data_count;

  die "# n col ids: ", scalar @col_ids, "; n dosages in row: ", $n_accessions, "\n" if($n_accessions != scalar @col_ids);
  my $marker_missing_data_fraction = $marker_missing_data_count/$n_accessions;

  if ($avg_pref_gt_prob >= $min_marker_avg_pref_gt_prob  and
      $marker_missing_data_fraction <= $max_marker_missing_data_fraction and
      $marker_minor_allele_frequency >= $min_marker_maf) { # there is good data for this marker; store and output later.
   # print STDERR "    maf: $marker_minor_allele_frequency  $min_marker_maf \n";
    push @row_ids, $row_id;	# store row (marker) id
    $rowid{$row_id} = 1; # also store in a hash to see whether all ids are distinct.
    push @alt_row_ids, $alt_row_id;
    $altrowid{$alt_row_id} = 1;
    for my $i (0..$ploidy) {
      # my $c = $marker_dosage_distribution[$i];
      $dosage_distribution[$i] +=  $marker_dosage_distribution[$i] // 0;
      #print STDERR $dosage_distribution[$i], " "; # /$n_elems_out, "  ";
    }
    $dosage_distribution[$ploidy+1] += $marker_missing_data_count;
    #print STDERR $dosage_distribution[$ploidy+1], ",  $marker_missing_data_count \n";

     while( my ($acc_idx, $the_dosage) = each @dosages_this_marker) {
	$accession_gt_strings[$acc_idx] .= "$the_dosage ";
    }
  }

  $markers_read_count++;
  print STDERR "# $markers_read_count lines of marker data read. \n" if($markers_read_count % 200 == 0);
}				# end loop over rows
close $fhin;
$info_string .= "# Done processing all $markers_read_count rows (markers).\n";
$info_string .= "# ", scalar @row_ids . " rows stored to be output.\n";
$info_string .= "# Data format: count.  ";
while (my($f, $c) = each %fieldused_count) {
  $info_string .= "$f: $c;  ";
}
$info_string .= "\n";
print $fhout $info_string;
print STDERR $info_string;

# #####  output  #####
# Use @rowids as marker ids, if they are all distinct,
if (scalar keys %rowid != scalar @row_ids) { # @row_ids are not all distinct; try @alt_row_ids
  if (scalar keys %altrowid == scalar @alt_row_ids) {
    @row_ids = @alt_row_ids;	# use @alt_row_ids 
  } else {	 # problem getting a set of unique marker identifiers.
    print STDERR scalar @row_ids, " ", scalar keys %rowid, " ", scalar @alt_row_ids, " ", scalar keys %altrowid, "\n";
    die "Neither row_ids nor alt_row_ids has all distinct identifiers.\n";
  }
}
my @row_ids_out = ();
for my $rid (@row_ids) {    # select the ids of the rows we are using.
  push @row_ids_out, $rid if(defined $rid);
}

my ($n_rows_out, $n_cols_out) = (scalar @row_ids_out, scalar @col_ids);
my $n_elements = $n_rows_out * $n_cols_out;
$info_string = '';
$info_string .= sprintf("# ploidy appears to be: %d\n", $ploidy);
$info_string .= sprintf("# outputting data for %d accessions and %d markers.\n", $n_cols_out, $n_rows_out);
$info_string .= sprintf("# dosage distribution:\n");
for my $i (0..$ploidy) {
  $info_string .= sprintf("#    %2d        %8d   %8.6f\n", $i, $dosage_distribution[$i], $dosage_distribution[$i]/$n_elements);
}
$info_string .= sprintf("# missing data %8d   %8.6f\n", $dosage_distribution[$ploidy+1], $dosage_distribution[$ploidy+1]/$n_elements);
print STDERR $info_string;
print $fhout $info_string; # if(! $plink_format);

# #####  print the output dosages  #####
# transposed
  print $fhout "MARKER ", join(" ", @row_ids_out), "\n";

while (my ($i, $col_id) = each @col_ids) {
    print $fhout "$col_id ";
    my $gtstr = $accession_gt_strings[$i];
    $gtstr =~ s/\t\s*$//;
    print $fhout "$gtstr\n";

}
close $fhout;
