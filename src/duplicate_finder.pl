#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw(min max sum);
use File::Spec 'splitpath';

# vcf file should have one of :
# dosage DS,
# genotype GT, e.g. 0/1 , 0/0/0/1  or  0|0|0|1
# allele depth AD, e.g. 21,19
# and if it has GP (genotype prob., e.g. 0.9,0.16,0.04) or GQ (genotype quality, e.g. 98) 
# we can also reject entries (i.e. regard as missing data) if these are not good enough

# usage:  vcf2dosage.pl -i <input vcf file>  -o <output file>


#my $transpose = 1; # default is to transpose; use -notrans to output untransposed.
# (duplicatesearch, find_parents require transposed)
# vcf: columns correspond to accessions, rows to markers

# if we don't believe can reliably resolve various heterozygous genotypes in polyploid case
# we can just lump together all heterozygous genotypes, map to just 3 genotypes:
my $map_to_012 = 0; # dosage = ploidy -> 2, 0 < dosage < ploidy -> 1, 0 -> 0, X -> X
my $field_to_use = 'GT'; # default GT. If requested field is not present exit (or ???)
# recognized choices are  GT (genotype e.g. '/1/0' ), DS (alternative allele dosage e.g. 2), AD (allele depths, e.g.'136:25' ).
### could add another option: GP, i.e. choose whichever gt has est. greatest probability.
my $ploidy = 2; # need to specify if ploidy > 2, will die if finds dosages greater than specified ploidy.-
my $inferred_ploidy = -1; # infer from data; die if > specified $ploidy
my $delta = 0.1; # if not $map_to_012, round to integer if within +- $delta
                 # if map_to_012 [0, $delta ->0], [1-$delta, $ploidy-1+$delta] -> 1, [$ploidy-$delta, $ploidy] -> 2
my $min_read_depth = 1;
my $vcf_filename = undef;
my $genotypes_filename = undef; # default: construct from input filename
my $missing_data_string = 'X';
my $minGQ = 0;			# if GQ present, must be >= this.
my $minGP = 0.0; # if GP present, there must be 1 genotype with prob >= $minGP; i.e. one genotype must be strongly preferred.
my $min_marker_avg_pref_gt_prob = -1.0; # default is negative (meaning don't filter on this)
my $max_marker_missing_data_fraction = 1.0; # remove markers with excessive missing data. Default is keep all.
my $min_marker_maf = 0;
my $max_distance = 0.25; 
my $info_string = "# command: " . join(" ", @ARGV) . "\n";
my $plink = 0;
my $use_alt_marker_ids = 0;

my $cluster_distance = 'auto'; # clusterer will attempt to choose a reasonable value.

GetOptions(
	   'input_file|vcf=s' => \$vcf_filename,
	   'output_file=s' => \$genotypes_filename,
	   #	   'GQmin=f' => \$minGQ,	# min genotype quality.
	   'GPmin=f' => \$minGP, # 
	   #	   'field=s' => \$field_to_use, # not implemented - just uses GT, so GT must be present in vcf file!
	   #	   'delta=f' => \$delta, 
	   #	   'min_read_depth=f' => \$min_read_depth,
	   'cluster_distance=f' => \$cluster_distance,
	   'dmax=f' => \$max_distance, 

	   'max_marker_md_fraction=f' => \$max_marker_missing_data_fraction,
	   'min_maf|maf_min=f' => \$min_marker_maf,
	   'plink!' => \$plink,
	   'alt_marker_ids!' => \$use_alt_marker_ids,

	   'ploidy=f' => \$ploidy,
	   'map_to_012!' => \$map_to_012,
	  );


my $clusterer_input_filename;
if (!defined $genotypes_filename) {
  (my $vol, my $dir, $genotypes_filename) = File::Spec->splitpath($vcf_filename);
  if ($genotypes_filename =~ /vcf$/) {
    $genotypes_filename =~ s/vcf$//; # remove vcf if present
    print STDERR "A: $genotypes_filename\n";
    $genotypes_filename =~ s/[.]$//; # remove final . if present
    print STDERR "B: $genotypes_filename\n";
  }
  $genotypes_filename .= "_gts";
}
print STDERR "genotypes_filename: $genotypes_filename \n";
my $vcf2gts_command = "vcf_to_gts -i $vcf_filename -p $minGP "; # for now uses GT field
$vcf2gts_command .= " -a " if($use_alt_marker_ids);
# "vcf2gts  -in $vcf_filename  -GQ $minGQ  -GP $minGP  -field $field_to_use ";


if ($plink) {	  #            *** analyze using plink ***
  $vcf2gts_command .= " -o $genotypes_filename ";
  print STDERR "vcf2gts command: $vcf2gts_command \n";
  system "$vcf2gts_command  -k";

  # get distances using plink
  my $plink_out_filename = $genotypes_filename . "_bin";
  print STDERR "plink_out_filename: $plink_out_filename \n"; #exit(0);
  my $plink_command1 = "plink1.9 --file $genotypes_filename --out $plink_out_filename --double-id --allow-extra-chr "; # --maf $min_marker_maf ";
  $plink_command1 .= " --maf $min_marker_maf " if($min_marker_maf > 0);
  system "$plink_command1"; # produces 3 files ending in .bed , .bin , and .fam
  my $plink_command2 = "plink1.9  --bfile $plink_out_filename --out $plink_out_filename --double-id --allow-extra-chr --distance-matrix ";
  $plink_command2 .= " --maf $min_marker_maf " if($min_marker_maf > 0);
  system "$plink_command2"; # produces files with endings .mdist (distance matrix), and .mdist.id (marker ids)
  my $cluster_filename_in = $plink_out_filename . ".dists";
  # filter out large distances and put in id1 id2 distance format 
  system "plnkout2dsout $plink_out_filename $cluster_filename_in $max_distance ";

  my $cluster_filename_out = $genotypes_filename . "_clusters";
  my $cluster_command = "clusterer -in $cluster_filename_in -out $cluster_filename_out -dcolumn 3 ";
  system "$cluster_command";

} else {       #                *** analyze using duplicate_search ***
  $vcf2gts_command .= " -o $genotypes_filename ";
  print STDERR "vcf_to_gts command: $vcf2gts_command \n";
  system "$vcf2gts_command";

  my $ds_distances_filename = $genotypes_filename . "_ds";
  my $ds_command = "duplicatesearch -i $genotypes_filename -a $min_marker_maf -e $max_distance -o $ds_distances_filename";
  print STDERR "duplicatesearch command: $ds_command\n";
  system "$ds_command";

  my $cluster_filename = $genotypes_filename . "_clusters";
  my $cluster_command = "clusterer -in $ds_distances_filename  -out $cluster_filename  -cluster_d $cluster_distance ";
  print STDERR "clusterer command: $cluster_command\n";
  system "$cluster_command";
}



