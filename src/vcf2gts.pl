#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw(min max sum);
use File::Spec 'splitpath';
use Devel::Size qw(size total_size);

# vcf file should have one of :

# genotype GT, e.g. 0/1  or  0|1
# dosage DS,
# allele depth AD, e.g. 21,19
# and if it has GP (genotype prob., e.g. 0.9,0.16,0.04) or GQ (genotype quality, e.g. 98) 
# we can also reject entries (i.e. regard as missing data) if these are not good enough

# usage:  vcf2gts.pl -i <input vcf file>  -o <output file>

# vcf: columns correspond to accessions, rows to markers

# if we don't believe can reliably resolve various heterozygous genotypes in polyploid case
# we can just lump together all heterozygous genotypes, map to just 3 genotypes:
my $map_to_012 = 0; # dosage = ploidy -> 2, 0 < dosage < ploidy -> 1, 0 -> 0, X -> X
my $field_to_use = 'GT'; # default GT. If requested field is not present exit (or ???)
# recognized choices are  GT (genotype e.g. '/1/0' ), DS (alternative allele dosage e.g. 2), AD (allele depths, e.g.'136:25' ).
my $ploidy = 2; # need to specify, will die if find dosages greater than specified ploidy.-
my $inferred_ploidy = -1; # infer from data; die if > specified $ploidy
my $delta = 0.1; # if not $map_to_012, round to integer if within +- $delta
                 # if map_to_012 [0, $delta ->0], [1-$delta, $ploidy-1+$delta] -> 1, [$ploidy-$delta, $ploidy] -> 2
my $min_read_depth = 5;
my $input_vcf_filename = undef;
my $output_genotypes_filename = undef; # default: construct from input filename
my $missing_data_string = 'X';
my $minGQ = 0;	 # if GQ present, must be >= this.
my $minGP = 0.0; # if GP present, there must be 1 genotype with prob >= $minGP; i.e. one genotype must be strongly preferred.
my $min_marker_avg_pref_gt_prob = -1.0; # default is negative (meaning don't filter on this)
my $max_marker_missing_data_fraction = 1.0; # remove markers with excessive missing data. Default is keep all.
my $min_marker_maf = 0;
my $info_string = "# command: " . join(" ", @ARGV) . "\n";
my $plink_format = 0;

GetOptions(
	   'input_file|vcf=s' => \$input_vcf_filename,
	   'output_file|dosage_file=s' => \$output_genotypes_filename,
	   'GQmin=f' => \$minGQ,	# min genotype quality.
	   'GPmin=f' => \$minGP,	# 
	   'field=s' => \$field_to_use,
	   'delta=f' => \$delta, 
	   'min_read_depth=f' => \$min_read_depth,
	   'max_marker_md_fraction=f' => \$max_marker_missing_data_fraction,
	   'min_maf=f' => \$min_marker_maf,
	   'ploidy=f' => \$ploidy,
	   'map_to_012!' => \$map_to_012,
	   'plink!' => \$plink_format,
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
if ($plink_format) {
  $output_genotypes_filename =~ s/dosages$//;
  $output_genotypes_filename =~ s/[.]$//;
  $output_markerids_filename = $output_genotypes_filename . ".map";
  $output_genotypes_filename .= ".ped";
  print STDERR "$output_markerids_filename\n";
}


open my $fhout, ">", "$output_genotypes_filename" or die "Couldn't open $output_genotypes_filename for writing.\n";
# ##########################################################

# #####  output run parameters  #####

$info_string .= "# output file: $output_genotypes_filename\n";
$info_string .= "# data field to use: $field_to_use\n";
$info_string .= "# delta: $delta ; min read depth: $min_read_depth\n";
$info_string .= "# GPmin: $minGP ; GQmin: $minGQ\n";
$info_string .= "# max marker missing data fraction: $max_marker_missing_data_fraction \n";
# $info_string .= "\n";
print $fhout $info_string if(! $plink_format);;
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
  my @cols = split(" ", $_);
  next if(/^\s*##/);
  if (/^\s*#/) {
    @col_ids = split(" ", $_);
    @col_ids = @col_ids[9..$#col_ids];
    last;
  }
}
print $fhout "# number of accession ids: #  ", scalar @col_ids, "\n" if(! $plink_format);
print STDERR "# number of accession ids: #  ", scalar @col_ids, "\n";

# #####  done reading accession ids  #####


# #####  read the rows with genotype data  #####
# print STDERR "plink map filename: $output_markerids_filename \n"; sleep(2);
open my $fhmap, ">", "$output_markerids_filename" or die "Couldn't open $output_markerids_filename for writing.\n" if($plink_format);
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
#  my $row_gts_str = '';
  my @cols = split(" ", $_);
  my $row_id = $cols[2]; # store row id in file (if these are not distinct, we will use $alt_row_id).
  my $alt_row_id = $cols[0] . "_" . $cols[1]; # construct a marker id from col[0] (chromosome number) and col[1] (position)
  my $chr_number = $cols[0];
  $chr_number =~ s/^chr//;
  my $plink_mapfile_str = "$chr_number\t$alt_row_id\t0\t" . $cols[1];
  #  print STDERR "$plink_mapfile_str\n";
  print $fhmap "$plink_mapfile_str\n" if($plink_format);
  my $format_str = $cols[8]; # e.g. 'GT:DS:AD' tells which types of data are present.
  my ($allele_zero, $allele_one) = @cols[3,4]; # e.g. 'A' and 'T' 

  my ($avg_pref_gt_prob, $pgtprob_count) = (0, 0); 

  # record which data fields are present (DS, GT, AD, GP, GQ)
  my @fields = split(':', $format_str);
  my $nfields = scalar @fields;
  my ($DSidx, $GTidx, $ADidx, $GPidx, $GQidx) = fields_present($format_str);

  # for this row (i.e. marker) loop over data for all accessions. Typical elements:
  # (GT:AD:DP): 0/0/0/0/0/1:79,18:97   ((hexaploid)genotype:allele depths:read depth)
  # (GT:DS:GP): 0|0:0.001:0.999,0.001,0    ((diploid)genotype:dosage:est_genotype_prob)
  @cols = @cols[9..$#cols];	# cols 9 and above have genotype info.
  # for my $e (@cols) {
  $n_accessions = 0;
  while (my($acc_idx, $e) = each @cols) {
    my $dosage = $missing_data_string; # indicates missing data
    my $first_allele = undef;	       # 0: ref allele, 1: alt allele.
    #   my $plnk_gt = "0\t0\t"; # this is what plink wants
    my @field_values = split(":", $e);
    die if(scalar @field_values != $nfields);

    # #####  if the data field is specified (GT, DS, GT, AD)  #####
    if ($field_to_use eq 'GT') {
      if ($GTidx >= 0) {	# GT; use genotype if present
	$fieldused_count{'GT'}++;
	if ($field_values[$GTidx] =~ /[.]/) { # no valid GT -> missing data
	} else {			      # count 0's and 1's
	  my $ref_allele_count = $field_values[$GTidx] =~ tr/0/0/;
	  my $alt_allele_count = $field_values[$GTidx] =~ tr/1/1/;
	  $first_allele = ($field_values[$GTidx] =~ /^(\S)/)? $1 : undef; #should be 0 or 1
	  $dosage = $alt_allele_count;
	}
      } else {
	die "Selected data field 'GT' not present.\n";
      }
    } elsif ($field_to_use eq 'DS') {
      if ($DSidx >= 0) {      # DS; use dosage if present in vcf file.
	$fieldused_count{'DS'}++;
	my $float_dosage = $field_values[$DSidx];
	my $int_dosage = int($float_dosage + 0.5); # round to integer
	if (abs($float_dosage - $int_dosage) <= $delta) { # float_dosage is close to integer, use
	  $dosage = $int_dosage;

	}			# else regard as missing data.
      } else {
	die "Selected data field 'DS' not present.\n";
      }
    } elsif ($field_to_use eq 'AD') {
      if ($ADidx >= 0) {
	# AD; use allele depth
	$fieldused_count{'AD'}++;
	my ($ref_depth, $alt_depth) = split(',', $field_values[$ADidx]);
	my $read_depth = $ref_depth + $alt_depth;
	if ($read_depth >= $min_read_depth) {
	  my $float_dosage = $ploidy*$alt_depth/$read_depth;
	  $dosage = int($float_dosage + 0.5);
	  if ($map_to_012) {
	    if (
		($float_dosage > $delta  and  $float_dosage < 1-$delta)
		or
		($float_dosage > $ploidy-1+$delta  and  $float_dosage < $ploidy-$delta)) {
	      $dosage = $missing_data_string;
	    }
	  } else {
	    $dosage = $missing_data_string if(
					      abs($float_dosage - $dosage) > $delta
					      or
					      ($float_dosage > $delta  and  $float_dosage < 1-$delta)
					      or
					      ($float_dosage > $ploidy-1+$delta  and  $float_dosage < $ploidy-$delta)
					     );
	  }
	}
      } else {
	die "Selected data field 'AD' not present.\n";
      }
    }

    ###   do quality checks here if GP or GQ available   ###
    if ($GPidx >= 0) {
      my @gt_probs = split(',', $field_values[$GPidx]);
      my $preferred_gt_prob = max(@gt_probs);
      # $avg_pref_gt_prob += $preferred_gt_prob;
      # $pgtprob_count++;
      $dosage = $missing_data_string if($preferred_gt_prob < $minGP);
    }
    if ($GQidx >= 0  and  $field_values[$GQidx] < $minGQ) {
      $dosage = $missing_data_string;
    }

    if ($dosage eq $missing_data_string) {
      $marker_missing_data_count++;
    } else {			# $dosage ne $missing_data_string
      $ploidy = $dosage if($dosage > $ploidy);
      $marker_dosage_distribution[$dosage]++;
      $marker_alt_allele_count += $dosage;
    }

    if ($plink_format) {
      my $plnk_gt = ($dosage eq $missing_data_string)? "0\t0\t" :
	($dosage == 0)? "$allele_zero\t$allele_zero\t" :
	($dosage == 1)? (($first_allele eq '0')? "$allele_zero\t$allele_one\t" : "$allele_one\t$allele_zero\t") :
	($dosage == 2)? "$allele_one\t$allele_one\t" : "0\t0\t";
      $accession_gt_strings[$acc_idx] .= $plnk_gt;
    } else {	      # not plink
      $accession_gt_strings[$acc_idx] .= "$dosage ";
    }
    $n_accessions++;
  }				# end loop over entries in a row
  die "Inferred ploidy ($inferred_ploidy) is greater than specified ploidy ($ploidy).\n" if($inferred_ploidy > $ploidy);
  my $marker_total_allele_count = $ploidy*($n_accessions - $marker_missing_data_count);
  my $marker_alt_allele_frequency = $marker_alt_allele_count/$marker_total_allele_count;
  my $marker_minor_allele_frequency = ($marker_alt_allele_frequency <= 0.5)?
    $marker_alt_allele_frequency :
    1.0 - $marker_alt_allele_frequency;

  $missing_data_count += $marker_missing_data_count;

  die "# n col ids: ", scalar @col_ids, "; n dosages in row: ", $n_accessions, "\n" if($n_accessions != scalar @col_ids);
  my $marker_missing_data_fraction = $marker_missing_data_count/$n_accessions;

  #  print STDERR "# avg pref gt prob: $avg_pref_gt_prob ;  marker md fraction: $marker_missing_data_fraction  $marker_missing_data_count\n";
  if ($avg_pref_gt_prob >= $min_marker_avg_pref_gt_prob  and
      $marker_missing_data_fraction <= $max_marker_missing_data_fraction and
      $marker_minor_allele_frequency >= $min_marker_maf) { # there is good data for this marker; store and output later.
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
  }

  $markers_read_count++;
  print STDERR "# $markers_read_count lines of marker data read. \n" if($markers_read_count % 200 == 0);
}				# end loop over rows
close $fhin;
close $fhmap if($plink_format);
$info_string .= "# Done processing all $markers_read_count rows (markers).\n";
$info_string .= "# ", scalar @row_ids . " rows stored to be output.\n";
$info_string .= "# Data format: count.  ";
while (my($f, $c) = each %fieldused_count) {
  $info_string .= "$f: $c;  ";
}
$info_string .= "\n";
print $fhout $info_string if(! $plink_format);
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
print $fhout $info_string if(! $plink_format);

# #####  print the output dosages  #####
# transposed
if (! $plink_format) {
  print $fhout "MARKER ", join(" ", @row_ids_out), "\n";
}
while (my ($i, $col_id) = each @col_ids) {
  # print col i (gts of ith accession) as a row
  if ($plink_format) {
    print $fhout "$col_id\t$col_id\t";
    print $fhout "0\t0\t0\t0\t"; # vcftools always puts these 4 zeroes - not sure why or whether necessary
  } else {
    print $fhout "$col_id ";
  }
  if (1 or $plink_format) {
    my $gtstr = $accession_gt_strings[$i];
    $gtstr =~ s/\t\s*$//;
    print $fhout "$gtstr\n";
  }
}
close $fhout;


sub fields_present{
  my $format_str = shift;
  my @fields = split(':', $format_str);
  my $nfields = scalar @fields;
  my ($DSidx, $GTidx, $ADidx, $GPidx, $GQidx) = (-1, -1, -1, -1, -1);
  #  my @fields_present = ();	# e.g. ('DS', 'GT', 'GP') 
  while (my($i, $f) = each @fields) {
    if ($f eq 'GT') { # genotype. e.g. 0/1 (unphased) or 0|1 (phased), or 0/0/0/1 (tetraploid)
      $GTidx = $i;
      #   push @fields_present, 'GT'; 
    } elsif ($f eq 'GP') { # genotype probability. e.g. 0.002,0.998,0.001   
      $GPidx = $i;
      #    push @fields_present, 'GP'; 
    } elsif ($f eq 'GQ') {	# genotype quality. e.g. 
      $GQidx = $i;
      #    push @fields_present, 'GT'; 
    } elsif ($f eq 'DS') { # dosage. e.g. 0 or 1, but can be non-integer.
      $DSidx = $i;
      #   push @fields_present, 'DS'; 
    } elsif ($f eq 'AD') {	# allele depth. e.g. 142,31
      $ADidx = $i;
      #    push @fields_present, 'AD'; 
    } 
  }
  return ($DSidx, $GTidx, $ADidx, $GPidx, $GQidx);
}
