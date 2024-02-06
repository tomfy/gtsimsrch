#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw (min max sum);

# simulate accessions with genotypes and phenotypes
# and add duplicates of accessions in that set
# with a bit of noise added to the genotypes (but not at present to phenotypes)

# then run plink to get .bed .bim .fam files from  .ped .map
# run gcta64 to do gwas analysis
# make first line of .loco.mlma file (gwas output) a comment, and
# remove lines with 'nan' in them

# my $genotypes_file = undef;
# my $phenotypes_file = undef;
# my $alleles_file = undef;
# my $map_file = undef;
# my $n_unique = 100;
# my $n_duplicate = 50;

my $simspec_file = '';
# the simspec file is something like e.g.:
# 499 Null1  0.05 0.95 0.0  0
# 2 Q1  0.05 0.95  0.08  0
# 499 Null2  0.05 0.95 0.0  0
# 499 Null3  0.05 0.95 0.0  0
# 2 Q2  0.05 0.95  0.07  0
# 499 Null4  0.05 0.95 0.0  0
# So, e.g. top line means:
# 499 markers, with ID Null1 maf between 0.05 and 0.095
# the next col is the size of the effect on the phenotypes (in some units I don't really understand)
# so zero for these 'Null' markers,
# and the next line specifies 2 markers with effect size 0.08


my $n_unique_acc = 600; # number of unique accessions simulated
my $n_dupe_acc = 400; # number of duplicate accessions
my $error_rate = 0.01; # when an accession is duplicated, errors at this rate are added to the duplicates and the original

my $n_reps = 10; # number of replicates, each one produces data sets with and without duplicates, and a gwas analysis
# of each of those. In particular there will be files with names like:
# u600_x.mlma  and u600d400_x.mlma where x is an integer ranging from 1 to $n_reps.

my $serial_number = 1; # starting serial number

 GetOptions(
	    'simspec_file=s' => \$simspec_file,
	    'n_acc|n_uniq=s' => \$n_unique_acc,
	    'n_dupe=s' => \$n_dupe_acc,
	    'error_rate=f' => \$error_rate,

	    'serial_number=i' => \$serial_number,
	    'n_reps=i' => \$n_reps,
	   );

for(1..$n_reps){
my $ufile = 'u' . $n_unique_acc . '_' . $serial_number;
my $udfile = 'u' . $n_unique_acc . 'd' . $n_dupe_acc . '_' . $serial_number;

# simulate; output is in form of  plink 'binary fileset', i.e. .fam, .bim, and .bed files.
system "plink1.9  --simulate-qt $simspec_file  --simulate-n $n_unique_acc  --out sim_temp ";

# convert to plink 'text fileset' so I can modify.
system "plink1.9  --bfile sim_temp  --recode  --out sim_temp ";

# make the phenotype file from .fam file
system "sel_columns '1,2,6' < sim_temp.fam > sim_temp.phen ";

# make the alleles file
system "sel_columns '2,5,6' < sim_temp.bim  > alleles ";

# make the duplicates, errors added to both unique and duplicate accessions
my $ped_file = "sim_temp.ped";
system "perl ~/gtsimsrch/src/rand_dupe_no_dupe.pl  $ped_file sim_temp.phen  alleles  $ufile $udfile $n_unique_acc $n_dupe_acc $error_rate ";

my $umap_file = $ufile . ".map";
my $map_file = "sim_temp.map";
system "ln -s $map_file $umap_file";
my $ubin_file = $ufile . "_bin";
system "plink1.9  --file $ufile  --out $ubin_file ";
my $uphen_file = $ufile . ".phen";

system "gcta64  --mlma --bfile $ubin_file  --pheno $uphen_file --out $ufile  --thread-num 2";

my $udmap_file = $udfile . ".map";
system "ln -s $map_file $udmap_file";
my $udbin_file = $udfile . "_bin";
system "plink1.9  --file $udfile  --out $udbin_file ";
my $udphen_file = $udfile . ".phen";

system "gcta64  --mlma --bfile $udbin_file  --pheno $udphen_file --out $udfile  --thread-num 2";

unlink $umap_file;
unlink $udmap_file;

$serial_number++;
}

