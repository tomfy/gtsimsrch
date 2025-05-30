#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename 'dirname';
use List::Util qw(min max sum);
use File::Spec 'splitpath';
use Time::HiRes qw( clock_gettime CLOCK_MONOTONIC );

use Cwd 'abs_path';
my ( $bindir, $libdir );
BEGIN {     # this has to go in Begin block so happens at compile time
  $bindir =
    dirname( abs_path(__FILE__) ) ; # the directory containing this script
}

# read in a vcf file, and process using either duplicate_search (default) or plink ( -plink );
# first run vcf_to_dsgm to convert vcf to format required by duplicate_search or plink;
# then run either duplicate_search or plink;
# then run clusterer to find clusters of accessions with near-identical genotype sets.

# vcf file should have :
# genotype GT, e.g. 0/1 , 0|0. (Diploid only)
# and can filter on GP (genotype prob., e.g. 0.9,0.16,0.04) ( e.g.  -GP 0.9 )
# would be nice to filter on GQ (genotype quality, e.g. 98) but not implemented yet. 

# usage:  duplicate_finder.pl -vcf <input vcf file>  -out <output file>

# duplicate_finder.pl  calls:
# vcf_to_dsgm
# duplicate_search (default) or ( -plink ) plink, plnkout2dsout
# clusterer.pl

{
my $command_str = "# command:  $0  " . join(" ", @ARGV);
my $field_to_use = 'GT'; # Presently GT is only option, must be present in vcf file.
# unimplemented alternatives: DS (alternative allele dosage e.g. 2), AD (allele depths, e.g.'136:25' ).

my $input_filename = undef;
my $ref_filename = undef;
my $genotypes_filename;
my $ref_genotypes_filename;

my $minGP = 0.9; # if GP present, there must be 1 genotype with prob >= $minGP; i.e. one genotype must be strongly preferred.
my $use_alt_marker_ids = 0; # default is use marker ids in col 3 of vcf file. -alt to construct marker ids from cols 1 and 2.

# my $minGQ = 0;	 # if GQ present, must be >= this. Not implemented.
# my $delta = 0.1; # if DS present, must be within $delta of an integer. Not implemented

my $plink = 0;
my $plnk2ds_perl = 0;	# default	# 0 -> C, 1 -> perl
my $plink_default_max_distance = 0.175;

my $chunk_size = 6;		# relevant only to duplicate_search
my $rng_seed = -1; # default: duplicate_search will get seed from clock
my $max_distance = 'default'; # duplicate_search only calculates distance if quick estimated distance is <= $max_distance; 'default' use duplicate_search's default.
# plink calculates all distances, and then we output only those <= $max_distance.
my $max_marker_missing_data_fraction = 0.25; # remove markers with excessive missing data.
my $max_accession_missing_data_fraction = 0.2; # Accessions with > missing data than this are excluded from analysis.
my $min_marker_maf = 0.08; # this is a good value for the yam 941 accession set.
my $full_duplicate_search_output = 1;
my $full_cluster_out = 1;
my $input_format = 'vcf';
my $ref_format = 'vcf';
my $histogram_path = 'histogram'; # might be e.g. 'perl /home/tomfy/Histogram_project/bin/histogram.pl'
my $histogram_agmr0 = 0;
my $histogram_color = undef;
my $nthreads = undef; # use default number of threads of duplicate_search or plink


my $graphics = 'gnuplot';
my $binwidth = 0.002;
my $plot_title = "counts of agmrs of accession pairs"; 
my $plot_xlabel = "agmr";
my $plot_ylabel = "counts";
my $histogram_filename = undef;

print "# duplicate_finder command: " . join(" ", @ARGV) . "\n";
my $distances_filename;
my $filename_stem;

# $cluster_distance defines how close accessions must be to put in same cluster.
my $cluster_distance = 'auto'; # default is 'auto': clusterer will attempt to choose a reasonable value.

GetOptions(

	   'input_filename|vcf_filename=s' => \$input_filename,

	  # 'format=s' => \$input_format, # either 'vcf' (default) or, if anything else -> dosage.
	   'output_file=s' => \$filename_stem,
	   'ref_filename|reference_filename=s' => \$ref_filename,
	  # 'ref_format=s' => \$ref_format, # default is 'vcf', anything else -> dosage 

	   # used by vcf_to_dsgm:
	   'min_gp|gp_min=f' => \$minGP, # if input file is vcf, and GP field present in vcf file, will filter on with this min value
	   'alt_marker_ids!' => \$use_alt_marker_ids, # -alt_marker_ids  -> construct marker id from chromosome number (col 1) and position (col 2)
	   #	   'GQmin=f' => \$minGQ,      # min genotype quality. Not implemented.
	   #       'delta=f' => \$delta,      # if

	   # to choose duplicate_search or plink:
	   'plink!' => \$plink,

	   # used by duplicate_search/plink:
	   'chunk_size|k=i' => \$chunk_size, # (duplicate_search only)
	   'seed|rand=i' => \$rng_seed,	     # (duplicate_search only)
	   'dmax|max_distance|distance_max=f' => \$max_distance,
	   'max_marker_md_fraction|max_marker_missing_data_fraction=f' => \$max_marker_missing_data_fraction,
	   'max_accession_md_fraction|accession_max_md_fraction=f' => \$max_accession_missing_data_fraction,
	   'min_maf|maf_min=f' => \$min_marker_maf,
	   'full_duplicate_search_output!' => $full_duplicate_search_output,
	   'threads|nthreads=i' => \$nthreads,

	   # used by clusterer:
	   'cluster_distance=f' => \$cluster_distance,
	   'full_cluster_output!' => \$full_cluster_out, #

	   'histogram_path=s' => \$histogram_path,
	   'show_agmr0!' => \$histogram_agmr0,
	   'xlabel=s' => \$plot_xlabel,
	   'binwidth|bw=s' => \$binwidth,
	   'graphics=s' => \$graphics,
	   'color=s' => \$histogram_color,
	   'histogram_filename=s' => \$histogram_filename,
	  );

my $t0 = clock_gettime(CLOCK_MONOTONIC);
my ($t1, $t2, $t3, $t4, $t5);
$input_format = determine_input_format($input_filename);
die "Input file $input_filename has unknown format.\n" if($input_format eq 'UNKNOWN');
if(defined $ref_filename){
  $ref_format = determine_input_format($ref_filename);
  die "Reference file $ref_filename has unknown format.\n" if($input_format eq 'UNKNOWN');
}

# $max_distance = -1 if($max_distance eq 'auto');

my $clusterer_input_filename;
if (!defined $filename_stem) {
  (my $vol, my $dir, $filename_stem) = File::Spec->splitpath($input_filename);
  # if ($filename_stem =~ /vcf$/) {
  #   $filename_stem =~ s/vcf$//; # remove vcf if present
  #   $filename_stem =~ s/[.]$//; # remove final . if present
  # }
}
$genotypes_filename = $filename_stem . "_dsgm";
# print  "# genotypes_filename: $genotypes_filename \n";
print  "# distances <= $max_distance will be found using ", ($plink)? "plink\n" : "duplicate_search\n";


if ($plink) {			#####  PLINK  #####
  $max_distance = $plink_default_max_distance if($max_distance eq 'default');
  my $plink_out_filename = $genotypes_filename . "_bin";

  my $plink_command1 = "plink1.9 --vcf $input_filename --double-id --out $filename_stem --vcf-min-gp $minGP ";
  $plink_command1 .= "  --threads $nthreads " if(defined $nthreads);
  print  "# plink command 1: $plink_command1\n";
  system "$plink_command1"; # produces 3 files ending in .bed , .bin , and .fam
 # my $plink_command2 = "plink1.9 --bfile $filename_stem --out $filename_stem --distance-matrix ";
   my $plink_command2 = "plink1.9 --bfile $filename_stem --out $filename_stem --distance 1-ibs flat-missing square";
  $plink_command2 .= " --maf $min_marker_maf --geno $max_marker_missing_data_fraction --mind $max_accession_missing_data_fraction ";
  $plink_command2 .= "  --threads $nthreads " if(defined $nthreads);
  print  "# plink command 2: $plink_command2\n";
  system "$plink_command2"; # produces files with endings .mdist (distance matrix), and .mdist.id (marker ids)

  $distances_filename = $filename_stem . ".dists";
  if ($plnk2ds_perl) {
    my $plnk2ds_abs_path = $bindir . "/plnkout2dsout.pl";
    system "$plnk2ds_abs_path  $filename_stem  $distances_filename $max_distance ";
  } else {		       # use C program to convert to ds format
    my $plnk_to_ds_abs_path = $bindir. "/plnkout_to_dsout ";
    my $distance_matrix_filename = $filename_stem . ".mdist";
    my $id_filename = $distance_matrix_filename . ".id";
    my $output_filename = $filename_stem . ".dists";
    system "$plnk_to_ds_abs_path  -i $id_filename  -d $distance_matrix_filename  -o $output_filename  -m $max_distance ";
  }

} else {			#####  DUPLICATE_SEARCH  #####
  $t1 = clock_gettime(CLOCK_MONOTONIC);
  print STDERR "# time to determine file format, etc.:  ", $t1 - $t0, "\n";
  if ($input_format eq 'vcf') {
    my $vcf2gts_command = $bindir . "/vcf_to_dsgm -input $input_filename -prob_min $minGP "; # for now uses GT field, can filter on GP
    $vcf2gts_command .= " -alternate_marker_ids " if($use_alt_marker_ids);
    $vcf2gts_command .= " -output $genotypes_filename ";
    print  "# vcf_to_dsgm command: $vcf2gts_command \n";
    print  "######### running vcf_to_dsgm ##########\n";
    system "$vcf2gts_command";
    print  "#########   vcf_to_dsgm done  ##########\n\n";
    $t2 = clock_gettime(CLOCK_MONOTONIC);
    print STDERR "# time for vcf_to_dsgm: ", $t2 - $t1, "\n";
  } else { # format was not 'vcf', assume input file has dosage format directly readable by duplicate_search
    $genotypes_filename = $input_filename;
  }

  if(defined $ref_filename  and  $ref_format eq 'vcf'){
      my $vcf2gts_command = $bindir . "/vcf_to_dsgm -input $ref_filename -pmin $minGP "; # for now uses GT field, can filter on GP
    $vcf2gts_command .= " -alternate_marker_ids " if($use_alt_marker_ids);
    $vcf2gts_command .= " -output $ref_genotypes_filename ";
    print  "# vcf_to_dsgm command: $vcf2gts_command \n";
    print  "######### running vcf_to_dsgm ##########\n";
      system "$vcf2gts_command";
    print  "#########   vcf_to_dsgm done  ##########\n\n";
  }else{
    $ref_genotypes_filename = $ref_filename;
  }
  
  $distances_filename = $filename_stem . ".dists";
  my $ds_command = $bindir . "/duplicate_search -input $genotypes_filename -maf_min $min_marker_maf -output $distances_filename";
  print STDERR "$ds_command\n";
  if ($max_distance ne 'default') {
    $ds_command .= " -distance_max $max_distance ";
  }
  if($full_duplicate_search_output){
    $ds_command .= " -format 2 ";
  }
  $ds_command .= " -ref $ref_genotypes_filename " if(defined $ref_filename);
  $ds_command .= " -marker_max_missing_data $max_marker_missing_data_fraction ";
  $ds_command .= " -chunk_size $chunk_size -accession_max_missing_data $max_accession_missing_data_fraction ";
  $ds_command .= " -seed $rng_seed " if($rng_seed > 0);
  $ds_command .= " -threads $nthreads " if(defined $nthreads);
  print  "# duplicate_search command: $ds_command\n";
  print  "######### running duplicate_search ##########\n";
  system "$ds_command";
  print  "#########  duplicate_search done  ##########\n\n";
  $t3 = clock_gettime(CLOCK_MONOTONIC);
  print STDERR "# time for duplicate_search: ", $t3 - $t2, "\n";
}

print  "######### running clusterer ##########\n";
my $cluster_filename = $filename_stem . "_clusters";
my $cluster_command = $bindir . "/clusterer.pl -in $distances_filename -out $cluster_filename -dcolumn 5 -cluster_d $cluster_distance ";
print  "# clusterer command: $cluster_command\n";
if (! $full_cluster_out) {
  $cluster_command .= " -nofull ";
}
# system "$cluster_command";
my $cluster_stdout = `$cluster_command`;
print "clusterer output to stdout: ", $cluster_stdout, "\n";
my $vline_xpos = undef;
if($cluster_stdout =~ /Max link distance: (\S+)/){
  $vline_xpos = $1;
}
 $t4 = clock_gettime(CLOCK_MONOTONIC);
  print STDERR "# time for clusterer: ", $t4 - $t3, "\n";
print  "#########  clusterer done  ##########\n\n";

print "#########  histogramming distances  ##########\n\n";

if(!defined $histogram_filename){ $histogram_filename = $filename_stem . '_distances_histogram.png'; }
my $histogram_command =
  #($histogram_agmr0  and  $full_duplicate_search_output)?
  #"$histogram_path -data $distances_filename:3/8 " :
  "$histogram_path -data '$distances_filename" . ':5"' . $plot_title . '"' . "' ";
 #if(lc $graphics eq 'gd');
$histogram_command .= " -vline $vline_xpos " if(defined $vline_xpos);
$histogram_command .= " -xlabel $plot_xlabel " if(defined $plot_xlabel);
$histogram_command .= " -ylabel $plot_ylabel " if(defined $plot_ylabel);
$histogram_command .= " -color $histogram_color " if(defined $histogram_color);
print STDERR "histogram command $histogram_command\n";
if(lc $graphics eq 'gnuplot'){
  $histogram_command .= ' -graphics gnuplot ';
}elsif(lc $graphics eq 'gd'){
  $histogram_command .= ' -graphics gd ';
  $histogram_command .= " -output $histogram_filename -bw $binwidth "; #  -png -noscreen -nointeractive "
}else{
  die "Graphics option $graphics is unknown. Options are 'gnuplot' and 'gd' (not case sensitive).\n";
}
# if($max_distance ne 'default'){
#   my $hi = $max_distance + 0.01;
#   $histogram_command .= " -hi $hi ";
# }
print "about to run command: [$histogram_command]\n";
system "$histogram_command";
print "#########  done histogramming  ##########\n\n";

}

sub determine_input_format{
  my $input_filename = shift;
  my $format = 'UNKNOWN';
  open my $fhin, "<", "$input_filename" or die "Couldn't open $input_filename for reading.\n";
  my $line;
  while($line = <$fhin>){
    if($line =~ /^#/){
      if($line =~ /^#CHROM/){
	$format = 'vcf';
	close $fhin;
	return $format;
      }
    }else{
      last;
    }
  }

  if($line =~ /^MARKER/){
    $format = 'dosage';
  }
  close $fhin;
  return $format;
}
