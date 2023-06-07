#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw(min max sum);
use File::Spec 'splitpath';

# read in a vcf file, and process using either duplicatesearch (default) or plink ( -plink );
# first run vcf_to_gts to convert vcf to format required by duplicatesearch or plink;
# then run either duplicatesearch or plink;
# then run clusterer to find clusters of accessions with near-identical genotype sets.

# vcf file should have :
# genotype GT, e.g. 0/1 , 0|0. (Diploid only)
# and can filter on GP (genotype prob., e.g. 0.9,0.16,0.04) ( e.g.  -GP 0.9 )
# would be nice to filter on GQ (genotype quality, e.g. 98) but not implemented yet. 

# usage:  duplicate_finder.pl -in <input vcf file>  -out <output file>

# duplicate_finder.pl  calls:
# vcf_to_gts
# duplicatesearch (default) or ( -plink ) plink, plnkout2dsout
# clusterer

my $field_to_use = 'GT'; # Presently GT is only option, must be present in vcf file.
# unimplemented alternatives: DS (alternative allele dosage e.g. 2), AD (allele depths, e.g.'136:25' ).

my $vcf_filename = undef;
my $ref_filename = undef;
#my $genotypes_filename = undef; # default: construct from input filename

my $minGP = 0.0; # if GP present, there must be 1 genotype with prob >= $minGP; i.e. one genotype must be strongly preferred.
my $use_alt_marker_ids = 0; # default is use marker ids in col 3 of vcf file. -alt to construct marker ids from cols 1 and 2.
# my $minGQ = 0;	 # if GQ present, must be >= this. Not implemented.
# my $delta = 0.1; # if DS present, must be within $delta of an integer. Not implemented

my $plink = 0;

my $chunk_size = 6; # relevant only to duplicatesearch
my $rng_seed = -1; # default duplicatesearch will get seed from clock
my $max_distance = 0.15; # duplicatesearch only calculates distance if quick estimated distance is <= $max_distance
my $max_marker_missing_data_fraction = 1.0; # remove markers with excessive missing data. Default is keep all.
my $max_accession_missing_data_fraction = 0.5; # Accessions with > missing data than this are excluded from analysis.
my $min_marker_maf = 0.01;
# plink calculates all distances, and then we output only those <= $max_distance.
my $info_string = "# command: " . join(" ", @ARGV) . "\n";
my $filename_stem;

my $oldway = 0;

# $cluster_distance defines how close accessions must be to put in same cluster.
my $cluster_distance = 'auto'; # default is 'auto': clusterer will attempt to choose a reasonable value.

GetOptions(
	   'input_file|vcf=s' => \$vcf_filename,
	   'output_file=s' => \$filename_stem,
	   'ref_filename|reference_filename=s' => \$ref_filename,

	   # used by vcf_to_gts:
	   'GPmin=f' => \$minGP, 
	   'alt_marker_ids!' => \$use_alt_marker_ids,
	   #	   'GQmin=f' => \$minGQ,      # min genotype quality. Not implemented.
	   #       'delta=f' => \$delta,      # if

	   # to choose duplicatesearch or plink:
	   'plink!' => \$plink,

	   # used by duplicatesearch/plink:
	   'chunk_size|k=i' => \$chunk_size, # (duplicatesearch only)
	    'seed|rand=i' => \$rng_seed, # (duplicatesearch only)
	   'dmax=f' => \$max_distance,
	   'max_marker_md_fraction|max_marker_missing_data_fraction=f' => \$max_marker_missing_data_fraction,
	   'max_accession_md_fraction|max_accession_md_fraction=f' => \$max_accession_missing_data_fraction,
	   'min_maf|maf_min=f' => \$min_marker_maf,

	   # used by clusterer:
	    'cluster_distance=f' => \$cluster_distance,
	  );


my $clusterer_input_filename;
if (!defined $filename_stem) {
  (my $vol, my $dir, $filename_stem) = File::Spec->splitpath($vcf_filename);
  if ($filename_stem =~ /vcf$/) {
    $filename_stem =~ s/vcf$//; # remove vcf if present
    $filename_stem =~ s/[.]$//; # remove final . if present
  }
}
my $genotypes_filename = $filename_stem . "_gts";
# print STDERR "# genotypes_filename: $genotypes_filename \n";
print STDERR "# distances <= $max_distance will be found using ", ($plink)? "plink\n" : "duplicatesearch\n";
my $vcf2gts_command = "vcf_to_gts -i $vcf_filename -p $minGP "; # for now uses GT field, can filter on GP
$vcf2gts_command .= " -a " if($use_alt_marker_ids);

if ($plink) {		      #            *** analyze using plink ***
    my $plink_out_filename = $genotypes_filename . "_bin";
  if($oldway){
  $vcf2gts_command .= " -o $genotypes_filename -k ";
  print STDERR "vcf conversion command: $vcf2gts_command \n";
  system "$vcf2gts_command";

  # get distances using plink

  print STDERR "plink_out_filename: $plink_out_filename \n"; #exit(0);
  my $plink_command1 = "plink1.9 --file $genotypes_filename --out $plink_out_filename --double-id --allow-extra-chr --make-bed "; # --maf $min_marker_maf ";
  $plink_command1 .= " --maf $min_marker_maf " if($min_marker_maf > 0);
  $plink_command1 .= " --geno $max_marker_missing_data_fraction  --mind $max_accession_missing_data_fraction ";
  print STDERR "# plink command 1: $plink_command1\n";

  system "$plink_command1"; # produces 3 files ending in .bed , .bin , and .fam
  my $plink_command2 = "plink1.9  --bfile $plink_out_filename --out $plink_out_filename --double-id --allow-extra-chr --distance-matrix ";
  $plink_command2 .= " --maf $min_marker_maf " if($min_marker_maf > 0);
  print STDERR "# plink command 2: $plink_command2\n";

  system "$plink_command2"; # produces files with endings .mdist (distance matrix), and .mdist.id (marker ids)
  my $cluster_filename_in = $genotypes_filename . ".dists";
  # filter out large distances and put in id1 id2 distance format 
  system "plnkout2dsout $plink_out_filename $cluster_filename_in $max_distance ";

  my $cluster_filename_out = $genotypes_filename . "_clusters";
  my $cluster_command = "clusterer -in $cluster_filename_in -out $cluster_filename_out -dcolumn 3 -cluster_d $cluster_distance ";
  system "$cluster_command";
}else{ # new way, using plink to process vcf file
  my $plink_command1 = "plink1.9 --vcf $vcf_filename --double-id --out $filename_stem --vcf-min-gp $minGP ";
  print STDERR "# plink command 1: $plink_command1\n";
  system "$plink_command1"; # produces 3 files ending in .bed , .bin , and .fam
  my $plink_command2 = "plink1.9 --bfile $filename_stem --out $filename_stem --distance-matrix ";
  $plink_command2 .= " --maf $min_marker_maf --geno $max_marker_missing_data_fraction --mind $max_accession_missing_data_fraction ";

  print STDERR "# plink command 2: $plink_command2\n";
  system "$plink_command2"; # produces files with endings .mdist (distance matrix), and .mdist.id (marker ids)

  my $cluster_filename_in = $filename_stem . ".dists";
  system "plnkout2dsout $filename_stem  $cluster_filename_in $max_distance ";

  my $cluster_filename_out = $filename_stem . "_clusters";
  my $cluster_command = "clusterer -in $cluster_filename_in -out $cluster_filename_out -dcolumn 3 -cluster_d $cluster_distance ";
  system "$cluster_command";
}

} else {       #                *** analyze using duplicate_search ***
  $vcf2gts_command .= " -o $genotypes_filename ";
  print STDERR "# vcf_to_gts command: $vcf2gts_command \n";
  print STDERR "######### running vcf_to_gts ##########\n";
  system "$vcf2gts_command";
  print STDERR "#########   vcf_to_gts done  ##########\n\n";

  my $ds_distances_filename = $filename_stem . ".dists";
  my $ds_command = "duplicatesearch -i $genotypes_filename -a $min_marker_maf -e $max_distance -o $ds_distances_filename";
  $ds_command .= " -r $ref_filename " if(defined $ref_filename);
  $ds_command .= " -x $max_marker_missing_data_fraction ";
  $ds_command .= " -k $chunk_size -f $max_accession_missing_data_fraction ";
  $ds_command .= " -s $rng_seed " if($rng_seed > 0);
  print STDERR "# duplicatesearch command: $ds_command\n";
  print STDERR "######### running duplicatesearch ##########\n";
  system "$ds_command";
  print STDERR "#########  duplicatesearch done  ##########\n\n";

  print STDERR "######### running clusterer ##########\n";
  my $cluster_filename = $filename_stem . "_clusters";
  my $cluster_command = "clusterer -in $ds_distances_filename -out $cluster_filename -dcolumn 3 -cluster_d $cluster_distance ";
  print STDERR "# clusterer command: $cluster_command\n";
  system "$cluster_command";
  print STDERR "#########  clusterer done  ##########\n\n";
}
