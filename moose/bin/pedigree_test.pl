#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(min max sum);
use Getopt::Long;
use File::Basename 'dirname';
use Time::HiRes qw( gettimeofday );
use Cwd 'abs_path';
my ( $bindir, $libdir );
BEGIN {     # this has to go in Begin block so happens at compile time
  $bindir =
    dirname( abs_path(__FILE__) ) ; # the directory containing this script
  $libdir = $bindir . '/../lib';
  $libdir = abs_path($libdir);	# collapses the bin/../lib to just lib
}
use lib $libdir;

use PedigreesDAG;
use Cluster1d;

# *****  run C program to get agmr, hgmr, r for parent-offspring pairs according to pedigree file,
# *****  then cluster with perl,
# *****  then run C program to test alternative pedigrees (if requested with -A option),
# *****  and analyze with perl to get best pedigrees.

# *****  input files:  *****
my $gtfilename = undef;
my $pedigree_table_filename = undef;

# *****  parameters for processing dosages, eliminating low-quality markers
# my $delta = 0.05; # real-number genotypes (dosages) are rounded to nearest integer (0,1,2) if within +- $delta
my $max_bad_gt_fraction = 0.2; # exclude markers with greater than this proportion of missing data.

# *****  output filenames:
my $base_output_filename = 'out';
my $output_pedigrees = 0;
my $output_genotype_matrix = 0;
my $c_output_filename_no_alt = 'pedigrees';
my $c_output_filename_alt = 'pedigrees_and_alternatives';
my $best_pedigrees_filename = 'best_pedigrees';

# *****  clustering control parameters:
my $pow = 'log';

# *****  control of search for alternative pedigrees:
my $find_alternatives = 0;


GetOptions(
	   # input filenames:
	   'gtsfile|gtfile|genotypesfile=s' => \$gtfilename,
	   'pedigreein|pedigreefile|pedtable=s' => \$pedigree_table_filename,
	   # control dosages -> cleaned genotypes
	#   'delta=f' => \$delta,
	   'max_bad_gt_fraction=f' => \$max_bad_gt_fraction,
	   # output filenames
	   'baseout|basename=s' => \$base_output_filename,
	   'outpedigrees!' => \$output_pedigrees,
	   'outmatrix!' => \$output_genotype_matrix,
	   'out_no_alt|out_noalt=s' => \$c_output_filename_no_alt,
	   'out_alt=s' => \$c_output_filename_alt,
	   'out_best=s' => \$best_pedigrees_filename,
	   # clustering control
	   'pow=s' => \$pow,
	   # control search for alternative pedigrees
	   'alternatives=i' => \$find_alternatives, # 0: none, 1: only for 'bad' pedigrees, >=2: do for all pedigrees.
	  );
die "No genotypes matrix filename provided.\n" if(!defined $gtfilename);


# *****  Read in the pedigree table:  ********
my $t_start = gettimeofday();
print STDERR "# Reading pedigrees from file: $pedigree_table_filename\n";
my $pedigrees = PedigreesDAG->new({pedigree_filename => $pedigree_table_filename});
my ($acyclic, $cyclic_id_string) = $pedigrees->is_it_acyclic();
print STDERR "# PedigreesDAG object created.  Is it acyclic: $acyclic\n";
if($acyclic == 0){
    print STDERR "# Warning: Pedigree file implies cycles in directed parent-offspring graph.\n";
    print STDERR "#          Childless accessions with ancestral cyclicities: $cyclic_id_string \n";
}
my $n_ids = scalar keys %{$pedigrees->id_node()};
my $n_without_offspring = scalar keys %{$pedigrees->childless_ids()};
my $n_with_offspring = $n_ids - $n_without_offspring;
print STDERR "# Number of ids in pedigree file: $n_ids\n",
    "# Number with/without offspring: $n_with_offspring / $n_without_offspring\n";

if ($output_pedigrees) {
  my $pedigree_output_filename = $base_output_filename . '_pedigrees';
  open my $fhout, ">", "$pedigree_output_filename";
  print $fhout $pedigrees->as_string();
  close $fhout;
}
printf(STDERR "# Time to read pedigree file, create PedigreesDAG obj: %5.3f\n\n", gettimeofday() - $t_start);
# *****  Done reading pedigree table. *********


# *****  Read in the genotype matrix file. Determine whether has dosages or 0,1,2,3 genotypes.
my $input_type = 'unknown';
my $c_program_gt_input_option = '-g';
open my $fhgt, "<", "$gtfilename";
while (my $line = <$fhgt>){
  last if($line =~ /^\s*MARKER/);
}
my $line = <$fhgt>;
my @cols = split(" ", $line);
if (scalar @cols == 2) {	# genotypes (0, 1, 2, or 3);
  if ($cols[1] =~ /[456789]+/) {
    print STDERR "# Should be only digits 0123 in genotypes. Exiting.\n";
    exit;
  }
  #  $input_type = 'genotypes0123';
  $c_program_gt_input_option = '-g';
} elsif (scalar @cols > 2) {	# dosages
  #  $input_type = 'dosages';
  $c_program_gt_input_option = '-d';
} else {
  print STDERR "# Number of columns is ", scalar @cols, ". Should be >= 2. Exiting.\n";
  exit;
}
# *****  Done loading the genotypes data  *******

# *****  Run c program to get stats on pedigrees in pedigree table:  *****
my $command = "pedigree_test $c_program_gt_input_option $gtfilename -p $pedigree_table_filename ";
$command .= " -x $max_bad_gt_fraction  -o $c_output_filename_no_alt ";
print STDERR "# Testing pedigrees using genotypes\n";
print STDERR  "# command: $command \n";
system "$command";
print STDERR "# after running c program.\n";
# *****  Test pedigrees in pedigree table:  *******************
open my $fhin, "<", "$c_output_filename_no_alt" or die "Couldn't open $c_output_filename_no_alt for reading.\n";
my @lines = <$fhin>;
print "# Number of pedigrees to be analyzed: ", scalar @lines, "\n";

# *****  Cluster agmr between parents in pedigree table
my @matpat_agmrs = ();
#my @array_of_lines_as_cols = ();
while (my ($j, $line) = each @lines) {
  my @cols = split(" ", $line);
  push @matpat_agmrs, $cols[5];
}
my $cluster1d_obj = Cluster1d->new({label => 'agmr between parents', xs => \@matpat_agmrs, pow => $pow});
my ($n_pts, $km_n_L, $km_n_R, $km_h_opt, $q, $kde_n_L, $kde_n_R, $kde_h_opt) = $cluster1d_obj->one_d_2cluster();
printf("# clustering of agmr between parents: %5d  k-means: %5d below %5d above %8.6f, q: %6.4f;  kde: %5d below %5d above %8.6f.\n", $n_pts, $km_n_L, $km_n_R, $km_h_opt, $q, $kde_n_L, $kde_n_R, $kde_h_opt);

my $max_self_agmr = $km_h_opt;

# #####  Clustering of hgmr, r, d1   ##########################
my @hgmr_denoms = ();
my @hgmrs = ();
my @r_denoms = ();
my @rs = ();
my @d_denoms = ();
my @ds = ();

while (my ($j, $line) = each @lines) {
  my @cols = split(" ", $line);

  my ($accid, $bad_gt_count) = @cols[0,1];
  my ($mat_id, $pat_id) = @cols[2,3];

  my $matpat_agmr = $cols[5];
  my ($mat_hgmr, $pat_hgmr) = @cols[7,11];
  my ($mat_r, $pat_r) = @cols[9,13];
#  my ($d1, $d2) = @cols[15,17];
  my $d = $cols[15];
  push @hgmr_denoms, @cols[6,10];
  push @hgmrs, ($mat_hgmr, $pat_hgmr);
  push @r_denoms, @cols[8,12];
  push @rs, ($mat_r, $pat_r);
  push @d_denoms, $cols[14];
  push @ds, $d;
}
my $median_matpat_agmr_denom = $matpat_agmrs[int(scalar @matpat_agmrs / 2)];
my $median_hgmr_denom = $hgmr_denoms[int(scalar @hgmr_denoms / 2)];
my $median_r_denom = $r_denoms[int(scalar @r_denoms / 2)];
my $median_d_denom = $d_denoms[int(scalar @d_denoms / 2)];

$cluster1d_obj = Cluster1d->new({label => 'hgmr', xs => \@hgmrs, pow => $pow});
($n_pts, $km_n_L, $km_n_R, $km_h_opt, $q, $kde_n_L, $kde_n_R, $kde_h_opt) = $cluster1d_obj->one_d_2cluster();
printf("# clustering of hgmr: %5d  k-means: %5d below %5d above %8.6f, q: %6.4f;  kde: %5d below %5d above %8.6f.\n", $n_pts, $km_n_L, $km_n_R, $km_h_opt, $q, $kde_n_L, $kde_n_R, $kde_h_opt);
my $max_ok_hgmr = $km_h_opt;

$cluster1d_obj = Cluster1d->new({label => 'r', xs => \@rs, pow => $pow});
($n_pts, $km_n_L, $km_n_R, $km_h_opt, $q, $kde_n_L, $kde_n_R, $kde_h_opt) = $cluster1d_obj->one_d_2cluster();
printf("# clustering of    r: %5d  k-means: %5d below %5d above %8.6f, q: %6.4f;  kde: %5d below %5d above %8.6f.\n", $n_pts, $km_n_L, $km_n_R, $km_h_opt, $q, $kde_n_L, $kde_n_R, $kde_h_opt);
my $max_self_r = $km_h_opt;

$cluster1d_obj = Cluster1d->new({label => 'd', xs => \@ds, pow => $pow});
($n_pts, $km_n_L, $km_n_R, $km_h_opt, $q, $kde_n_L, $kde_n_R, $kde_h_opt) = $cluster1d_obj->one_d_2cluster();
printf("# clustering of   d: %5d  k-means: %5d below %5d above %8.6f, q: %6.4f;  kde: %5d below %5d above %8.6f.\n", $n_pts, $km_n_L, $km_n_R, $km_h_opt, $q, $kde_n_L, $kde_n_R, $kde_h_opt);
my $max_ok_d = $km_h_opt;

# *****  Look at alternative pedigrees if requested
if ($find_alternatives > 0) {
  $max_self_agmr = sprintf("%.6f", $max_self_agmr);
  $max_ok_hgmr = sprintf("%.6f", $max_ok_hgmr);
  $max_self_r = sprintf("%.6f", $max_self_r);
  $max_ok_d = sprintf("%.6f", $max_ok_d);
    
#  $gtfilename = 'genotype_matrix_out';
  $command = "pedigree_test -g $gtfilename -p $pedigree_table_filename   -x $max_bad_gt_fraction -A $find_alternatives -o $c_output_filename_alt ";
  $command .= " -a $max_self_agmr  -h $max_ok_hgmr -r $max_self_r -D $max_ok_d";
  print "# command: $command \n";
  system "$command";


  open my $fh_alt, "<", "$c_output_filename_alt";
  @lines = <$fh_alt>;

  my $factor = 0.33;
  my %category_counts = ();
  my $bad_denoms_count = 0;

  open my $fhbest, ">", "$best_pedigrees_filename";
  while (my ($j, $line) = each @lines) {
    my @cols = split(" ", $line);
    my ($acc_id, $acc_md) = @cols[0,1];
    my ($mat_id, $pat_id) = @cols[2,3];
    my $id_pair = "$mat_id $pat_id";
    my $ped_d = $cols[15];
    my $category_string = '';
    my %allped_d = ();
    my $ok_pedigrees_count = 0; # counts all ok (small d) pedigrees for this accession, both pedigree from table and alternatives.
    my $denoms_ok = are_denoms_ok(\@cols, 4, $factor, $median_matpat_agmr_denom, $median_hgmr_denom, $median_r_denom, $median_d_denom, $factor);
    $category_string = ($denoms_ok)? category(\@cols, 4, $max_self_agmr, $max_ok_hgmr, $max_self_r, $max_ok_d) : 'x xx xx x';
    my $ped_str = sprintf("  ped  %20s %20s  ", $mat_id, $pat_id) . "  $category_string";
    if ($denoms_ok) {
      $ok_pedigrees_count++ if ($category_string eq '0 00 00 0'  or  $category_string eq '1 01 01 0');
      $category_counts{$category_string}++;
    } else {
      $bad_denoms_count++;
    }
    $allped_d{$ped_str} = $ped_d;

    # alternative pedigrees:
    my $n_alternatives = $cols[16] // 0; #
    my %okalt_d = ();
    for (my $i = 0; $i < $n_alternatives; $i++) {
      my $first = 19 + $i*14;
      my $alt_category_string = '';
      my $alt_id_pair = sprintf("%20s %20s", $cols[$first-2],  $cols[$first-1]);
      my $alt_d = $cols[$first+11];
      $denoms_ok = are_denoms_ok(\@cols, $first, $factor, $median_matpat_agmr_denom, $median_hgmr_denom, $median_r_denom, $median_d_denom, $factor);
      if ($denoms_ok) {
	$alt_category_string = category(\@cols, $first, $max_self_agmr, $max_ok_hgmr, $max_self_r, $max_ok_d);
	if ($alt_category_string eq '0 00 00 0'  or $alt_category_string eq '1 01 01 0') {
	  $ok_pedigrees_count++;
	  $okalt_d{"  alt  $alt_id_pair  $alt_category_string"} = $alt_d;
	}
      } else {
	$alt_category_string = 'x xx xx x';
      }
      $allped_d{"  alt  $alt_id_pair  $alt_category_string"} = $alt_d;
    }

    
    my $output_string = $cols[0] . "  $ok_pedigrees_count";
    my @sorted_alts = #sort {$okalt_d{$a} <=> $okalt_d{$b} } keys %okalt_d;
    sort { # sort first by number of x's (low to high), then by d (low_to_high)
  my $acat = ($a =~ /(\S+\s+\S+\s+\S+\s+\S+)\s*$/)? $1 : 'x xx xx x';
  my $bcat = ($b =~ /(\S+\s+\S+\s+\S+\s+\S+)\s*$/)? $1 : 'x xx xx x';
  my $ax = () = $acat =~ /x/g;
  my $bx = () = $bcat =~ /x/g;
  (($ax <=> $bx) or ($okalt_d{$a} <=> $okalt_d{$b})); } keys %okalt_d;
    
    my @sorted_allpeds = # sort {$allped_d{$a} <=> $allped_d{$b} } keys %allped_d;
 sort {
  my $acat = ($a =~ /(\S+\s+\S+\s+\S+\s+\S+)\s*$/)? $1 : 'x xx xx x';
  my $bcat = ($b =~ /(\S+\s+\S+\s+\S+\s+\S+)\s*$/)? $1 : 'x xx xx x';
  my $ax = () = $acat =~ /x/g;
  my $bx = () = $bcat =~ /x/g;
  (($ax <=> $bx) or ($allped_d{$a} <=> $allped_d{$b})); } keys %allped_d;
      
    if ($category_string eq '0 00 00 0'  or  $category_string eq '1 01 01 0') { # pedigree from table is 'good'
      $output_string .= "  ped" . $ped_str . "  " . $ped_d; # ped  " . $id_pair . "  $category_string  $ped_d";
      if (scalar @sorted_alts > 0) {
	$output_string .= $sorted_alts[0] . "  " . $okalt_d{$sorted_alts[0]};
      }
    } elsif (scalar @sorted_alts > 0) {
      $output_string .= "  alt";
      for (my $i = 0; $i < min(scalar @sorted_alts, 2); $i++) {
	$output_string .= $sorted_alts[$i] . "  " . $okalt_d{$sorted_alts[$i]};
      }
    } else {
      $output_string .= "  none";
      for (my $i = 0; $i < min(scalar @sorted_allpeds, 2); $i++) {
	$output_string .= $sorted_allpeds[$i] . "  " . $allped_d{$sorted_allpeds[$i]};
      }
    }
    $output_string .= "\n";
    print $fhbest $output_string; # output 1 line about likely pedigrees for this accession.
  }
  print $fhbest "# Bad denoms count: $bad_denoms_count\n";
  my @scategories = sort {$a cmp $b} keys %category_counts;
  # while (my($cat, $count) = each %category_counts) {
  for my $cat (@scategories) {
    print "$cat ", $category_counts{$cat}, "\n";
  }
}


sub are_denoms_ok{
  my @cols = @{my $cls = shift};
  my $first = shift;
  my $factor = shift;
  my ($median_matpat_agmr_denom, $median_hgmr_denom, $median_r_denom, $median_d_denom) = @_;
  my $last = $first + 11;
  my ($FMagmr_denom, $FMagmr, $Fhgmr_denom, $Fhgmr, $Fr_denom, $Fr, $Mhgmr_denom, $Mhgmr, $Mr_denom, $Mr, $d_denom, $d) = @cols[$first..$last];
  my $agmr_denom_ok = ($FMagmr_denom >= $factor*$median_matpat_agmr_denom);
  my $hgmr_denoms_ok = ($Fhgmr_denom >= $factor*$median_hgmr_denom  and  $Mhgmr_denom >= $factor*$median_hgmr_denom);
  my $r_denoms_ok = ($Fr_denom >= $factor*$median_r_denom  and  $Mr_denom >= $factor*$median_r_denom);
  my $d_denom_ok = ($d_denom >= $factor*$median_d_denom);
  my $denoms_ok = ($agmr_denom_ok  and  $hgmr_denoms_ok  and  $r_denoms_ok  and  $d_denom_ok);
  return $denoms_ok;
}

sub category{
  my @cols = @{my $cls = shift};
  my $first = shift;
  my ($max_self_agmr, $max_ok_hgmr, $max_self_r, $max_ok_d) = @_;
  my $last = $first + 11;
  my ($FMagmr_denom, $FMagmr, $Fhgmr_denom, $Fhgmr, $Fr_denom, $Fr, $Mhgmr_denom, $Mhgmr, $Mr_denom, $Mr, $d_denom, $d) = @cols[$first..$last];
  my $category_string = '';
  $category_string .= ($FMagmr <= $max_self_agmr)? '0' : '1';
  $category_string .= ' ';
  $category_string .= ($Fhgmr <= $max_ok_hgmr)? '0' : '1';
  $category_string .= ($Fr <= $max_self_r)? '0' : '1';
  $category_string .= ' ';
  $category_string .= ($Mhgmr <= $max_ok_hgmr)? '0' : '1';
  $category_string .= ($Mr <= $max_self_r)? '0' : '1';
  $category_string .= ' ';
  $category_string .= ($d <= $max_ok_d)? '0' : '1';
  return $category_string;
}

