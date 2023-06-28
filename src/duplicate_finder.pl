#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename 'dirname';
use List::Util qw(min max sum);
use File::Spec 'splitpath';

use Cwd 'abs_path';
my ( $bindir, $libdir );
BEGIN {     # this has to go in Begin block so happens at compile time
  $bindir =
    dirname( abs_path(__FILE__) ) ; # the directory containing this script
}

# read in a vcf file, and process using either duplicatesearch (default) or plink ( -plink );
# first run vcf_to_gts to convert vcf to format required by duplicatesearch or plink;
# then run either duplicatesearch or plink;
# then run clusterer to find clusters of accessions with near-identical genotype sets.

# vcf file should have :
# genotype GT, e.g. 0/1 , 0|0. (Diploid only)
# and can filter on GP (genotype prob., e.g. 0.9,0.16,0.04) ( e.g.  -GP 0.9 )
# would be nice to filter on GQ (genotype quality, e.g. 98) but not implemented yet. 

# usage:  duplicate_finder.pl -vcf <input vcf file>  -out <output file>

# duplicate_finder.pl  calls:
# vcf_to_gts
# duplicatesearch (default) or ( -plink ) plink, plnkout2dsout
# clusterer.pl

my $field_to_use = 'GT'; # Presently GT is only option, must be present in vcf file.
# unimplemented alternatives: DS (alternative allele dosage e.g. 2), AD (allele depths, e.g.'136:25' ).

my $vcf_filename = undef;
my $ref_filename = undef;

my $minGP = 0.9; # if GP present, there must be 1 genotype with prob >= $minGP; i.e. one genotype must be strongly preferred.
my $use_alt_marker_ids = 0; # default is use marker ids in col 3 of vcf file. -alt to construct marker ids from cols 1 and 2.

# my $minGQ = 0;	 # if GQ present, must be >= this. Not implemented.
# my $delta = 0.1; # if DS present, must be within $delta of an integer. Not implemented

my $plink = 0;
my $use_plnk2ds = 0;
my $plink_default_max_distance = 0.175;

my $chunk_size = 6; # relevant only to duplicatesearch
my $rng_seed = -1;  # default: duplicatesearch will get seed from clock
my $max_distance = 'auto'; # duplicatesearch only calculates distance if quick estimated distance is <= $max_distance; 'auto' -> get random sample of dists, use to choose $max_distance.
# plink calculates all distances, and then we output only those <= $max_distance.
my $max_marker_missing_data_fraction = 0.25; # remove markers with excessive missing data.
my $max_accession_missing_data_fraction = 0.5; # Accessions with > missing data than this are excluded from analysis.
my $min_marker_maf = 0.08; # this is a good value for the yam 941 accession set.
my $full_cluster_out = 1;

print "# duplicate_finder command: " . join(" ", @ARGV) . "\n";
my $distances_filename;
my $filename_stem;

# $cluster_distance defines how close accessions must be to put in same cluster.
my $cluster_distance = 'auto'; # default is 'auto': clusterer will attempt to choose a reasonable value.

GetOptions(
	   'vcf|input_filename=s' => \$vcf_filename,
	   'output_file=s' => \$filename_stem,
	   'ref_filename|reference_filename=s' => \$ref_filename,

	   # used by vcf_to_gts:
	   'min_gp|gp_min=f' => \$minGP, 
	   'alt_marker_ids!' => \$use_alt_marker_ids,
	   #	   'GQmin=f' => \$minGQ,      # min genotype quality. Not implemented.
	   #       'delta=f' => \$delta,      # if

	   # to choose duplicatesearch or plink:
	   'plink!' => \$plink,

	   # used by duplicatesearch/plink:
	   'chunk_size|k=i' => \$chunk_size, # (duplicatesearch only)
	   'seed|rand=i' => \$rng_seed,	     # (duplicatesearch only)
	   'dmax|max_distance=f' => \$max_distance,
	   'max_marker_md_fraction|max_marker_missing_data_fraction=f' => \$max_marker_missing_data_fraction,
	   'max_accession_md_fraction|accession_max_md_fraction=f' => \$max_accession_missing_data_fraction,
	   'min_maf|maf_min=f' => \$min_marker_maf,

	   # used by clusterer:
	   'cluster_distance=f' => \$cluster_distance,
	   'full_cluster_output!' => \$full_cluster_out, #
	  );

# $max_distance = -1 if($max_distance eq 'auto');

my $clusterer_input_filename;
if (!defined $filename_stem) {
  (my $vol, my $dir, $filename_stem) = File::Spec->splitpath($vcf_filename);
  if ($filename_stem =~ /vcf$/) {
    $filename_stem =~ s/vcf$//; # remove vcf if present
    $filename_stem =~ s/[.]$//; # remove final . if present
  }
}
my $genotypes_filename = $filename_stem . "_gts";
# print  "# genotypes_filename: $genotypes_filename \n";
print  "# distances <= $max_distance will be found using ", ($plink)? "plink\n" : "duplicatesearch\n";


if ($plink) {  #####  PLINK  #####
  $max_distance = $plink_default_max_distance if($max_distance eq 'auto');
  my $plink_out_filename = $genotypes_filename . "_bin";

  my $plink_command1 = "plink1.9 --vcf $vcf_filename --double-id --out $filename_stem --vcf-min-gp $minGP ";
  print  "# plink command 1: $plink_command1\n";
  system "$plink_command1"; # produces 3 files ending in .bed , .bin , and .fam
  my $plink_command2 = "plink1.9 --bfile $filename_stem --out $filename_stem --distance-matrix ";
  $plink_command2 .= " --maf $min_marker_maf --geno $max_marker_missing_data_fraction --mind $max_accession_missing_data_fraction ";

  system "$plink_command2"; # produces files with endings .mdist (distance matrix), and .mdist.id (marker ids)

  $distances_filename = $filename_stem . ".dists";
  if ($use_plnk2ds) {
    my $plnk2ds_abs_path = $bindir . "/plnkout2dsout.pl";
    system "$plnk2ds_abs_path  $filename_stem  $distances_filename $max_distance ";
  } else {		       # use C program to convert to ds format
    my $plnk_to_ds_abs_path = $bindir. "/plnkout_to_dsout ";
    my $distance_matrix_filename = $filename_stem . ".mdist";
    my $id_filename = $distance_matrix_filename . ".id";
    my $output_filename = $filename_stem . ".dists";
    system "$plnk_to_ds_abs_path  -i $id_filename  -d $distance_matrix_filename  -o $output_filename  -m $max_distance ";
  }

} else {   #####  DUPLICATESEARCH  #####

  my $vcf2gts_command = $bindir . "/vcf_to_gts -input $vcf_filename -pmin $minGP "; # for now uses GT field, can filter on GP
  $vcf2gts_command .= " -alternate_marker_ids " if($use_alt_marker_ids);
  $vcf2gts_command .= " -output $genotypes_filename ";
  print  "# vcf_to_gts command: $vcf2gts_command \n";
  print  "######### running vcf_to_gts ##########\n";
  system "$vcf2gts_command";
  print  "#########   vcf_to_gts done  ##########\n\n";

  $distances_filename = $filename_stem . ".dists";
  my $ds_command = $bindir . "/duplicatesearch -input $genotypes_filename -maf_min $min_marker_maf -output $distances_filename";
  if ($max_distance ne 'auto') {
    $ds_command .= " -max_est_distance $max_distance ";
  }
  $ds_command .= " -ref $ref_filename " if(defined $ref_filename);
  $ds_command .= " -marker_max_missing_data $max_marker_missing_data_fraction ";
  $ds_command .= " -chunk_size $chunk_size -accession_max_missing_data $max_accession_missing_data_fraction ";
  $ds_command .= " -seed $rng_seed " if($rng_seed > 0);
  print  "# duplicatesearch command: $ds_command\n";
  print  "######### running duplicatesearch ##########\n";
  system "$ds_command";
  print  "#########  duplicatesearch done  ##########\n\n";
}

 print  "######### running clusterer ##########\n";
  my $cluster_filename = $filename_stem . "_clusters";
  my $cluster_command = $bindir . "/clusterer.pl -in $distances_filename -out $cluster_filename -dcolumn 3 -cluster_d $cluster_distance ";
print  "# clusterer command: $cluster_command\n";
if(! $full_cluster_out){
  $cluster_command .= " -nofull ";
}
  system "$cluster_command";
  print  "#########  clusterer done  ##########\n\n";
