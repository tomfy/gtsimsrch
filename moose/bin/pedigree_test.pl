#!/usr/bin/perl
use strict;
use warnings;
use Scalar::Util qw(looks_like_number);
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
my $max_bad_gt_fraction = 0.1; # exclude markers with greater than this proportion of missing data.

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
# my $input_type = 'unknown';
# open my $fhgt, "<", "$gtfilename";
# while (my $line = <$fhgt>){
#   last if($line =~ /^\s*MARKER/);
# }
# my $line = <$fhgt>;
# my @cols = split(" ", $line);
# if (scalar @cols == 2) {	# genotypes (0, 1, 2, or 3);
#   if ($cols[1] =~ /[456789]+/) {
#     print STDERR "# Should be only digits 0123 in genotypes. Exiting.\n";
#     exit;
#   }
#   #  $input_type = 'genotypes0123';
# } elsif (scalar @cols > 2) {	# dosages
#   #  $input_type = 'dosages';
#   $c_program_gt_input_option = '-d';
# } else {
#   print STDERR "# Number of columns is ", scalar @cols, ". Should be >= 2. Exiting.\n";
#   exit;
# }
# *****  Done loading the genotypes data  *******

# *****  Run c program to get stats on pedigrees in pedigree table:  *****
my $command = "pedigree_test -g $gtfilename -p $pedigree_table_filename ";
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
  next if($line =~ /^\s*#/);
  my @cols = split(" ", $line);
  if(looks_like_number($cols[5])){
      push @matpat_agmrs, $cols[5];
  }
}
my $cluster1d_obj = Cluster1d->new({label => 'agmr between parents', xs => \@matpat_agmrs, pow => $pow}); #, median_denom => $median_matpat_agmr_denom});
my ($n_pts, $km_n_L, $km_n_R, $km_h_opt, $q, $kde_n_L, $kde_n_R, $kde_h_opt) = $cluster1d_obj->one_d_2cluster();
printf("# clustering of agmr between parents: %5d  k-means: %5d below %5d above %8.6f, q: %6.4f;  kde: %5d below %5d above %8.6f.\n", $n_pts, $km_n_L, $km_n_R, $km_h_opt, $q, $kde_n_L, $kde_n_R, $kde_h_opt);

my $max_self_agmr = $km_h_opt;

# exit;

# #####  Clustering of hgmr, r, z, d  ##########################
my @hgmr_denoms = ();
my @hgmrs = ();
my @r_denoms = ();
my @rs = ();
my @z_denoms = ();
my @zs = ();
my @d_denoms = ();
my @ds = ();

while (my ($j, $line) = each @lines) {
  next if($line =~ /^\s*#/);
  my @cols = split(" ", $line);

  my ($accid, $bad_gt_count) = @cols[0,1];
  my ($mat_id, $pat_id) = @cols[2,3];

#  my $matpat_agmr = $cols[5];
#  my ($mat_hgmr, $pat_hgmr) = @cols[7,11];
#  my ($mat_r, $pat_r) = @cols[9,13];
#  my ($d1, $d2) = @cols[15,17];
#  my $z = $cols[15];
#  my $d = $cols[17];

  my ($mat_hgmr_denom, $mat_hgmr) = @cols[6,7];
  if(looks_like_number($mat_hgmr)){
      push @hgmr_denoms, $mat_hgmr_denom;
      push @hgmrs, $mat_hgmr;
  }
   my ($pat_hgmr_denom, $pat_hgmr) = @cols[10,11];
  if(looks_like_number($pat_hgmr)){
      push @hgmr_denoms, $pat_hgmr_denom;
      push @hgmrs, $pat_hgmr;
  }

  my ($mat_R_denom, $mat_R) = @cols[8,9];
  if(looks_like_number($mat_R)){
      push @r_denoms, $mat_R_denom;
      push @rs, $mat_R;
  }
   my ($pat_R_denom, $pat_R) = @cols[12,13];
  if(looks_like_number($pat_R)){
      push @r_denoms, $pat_R_denom;
      push @rs, $pat_R;
  }

   my ($z_denom, $z) = @cols[14,15];
  if(looks_like_number($z)){
      push @z_denoms, $z_denom;
      push @zs, $z;
  }
  
   my ($d_denom, $d) = @cols[16,17];
  if(looks_like_number($d)){
      push @d_denoms, $d_denom;
      push @ds, $d;
  }

# push @hgmr_denoms, @cols[6,10];
#  push @hgmrs, ($mat_hgmr, $pat_hgmr);
#  push @r_denoms, @cols[8,12];
#  push @rs, ($mat_r, $pat_r);
#  push @z_denoms, $cols[14];
#  push @zs, $z;
#  push @d_denoms, $cols[16];
#  push @ds, $d;
}
my $median_matpat_agmr_denom = $matpat_agmrs[int(scalar @matpat_agmrs / 2)];
my $median_hgmr_denom = $hgmr_denoms[int(scalar @hgmr_denoms / 2)];
my $median_r_denom = $r_denoms[int(scalar @r_denoms / 2)];
my $median_z_denom = $z_denoms[int(scalar @z_denoms / 2)];
my $median_d_denom = $d_denoms[int(scalar @d_denoms / 2)];

$cluster1d_obj = Cluster1d->new({label => 'hgmr', xs => \@hgmrs, pow => $pow}); #
($n_pts, $km_n_L, $km_n_R, $km_h_opt, $q, $kde_n_L, $kde_n_R, $kde_h_opt) = $cluster1d_obj->one_d_2cluster();
printf("# clustering of hgmr: %5d  k-means: %5d below %5d above %8.6f, q: %6.4f;  kde: %5d below %5d above %8.6f.\n", $n_pts, $km_n_L, $km_n_R, $km_h_opt, $q, $kde_n_L, $kde_n_R, $kde_h_opt);
my $max_ok_hgmr = $km_h_opt;

$cluster1d_obj = Cluster1d->new({label => 'r', xs => \@rs, pow => $pow});
($n_pts, $km_n_L, $km_n_R, $km_h_opt, $q, $kde_n_L, $kde_n_R, $kde_h_opt) = $cluster1d_obj->one_d_2cluster();
printf("# clustering of    r: %5d  k-means: %5d below %5d above %8.6f, q: %6.4f;  kde: %5d below %5d above %8.6f.\n", $n_pts, $km_n_L, $km_n_R, $km_h_opt, $q, $kde_n_L, $kde_n_R, $kde_h_opt);
my $max_self_r = $km_h_opt;

$cluster1d_obj = Cluster1d->new({label => 'z', xs => \@zs, pow => $pow});
($n_pts, $km_n_L, $km_n_R, $km_h_opt, $q, $kde_n_L, $kde_n_R, $kde_h_opt) = $cluster1d_obj->one_d_2cluster();
printf("# clustering of   z: %5d  k-means: %5d below %5d above %8.6f, q: %6.4f;  kde: %5d below %5d above %8.6f.\n", $n_pts, $km_n_L, $km_n_R, $km_h_opt, $q, $kde_n_L, $kde_n_R, $kde_h_opt);
my $max_ok_z = $km_h_opt;

print STDERR "Clustering d!!!!!\n";
$cluster1d_obj = Cluster1d->new({label => 'd', xs => \@ds, pow => $pow});
($n_pts, $km_n_L, $km_n_R, $km_h_opt, $q, $kde_n_L, $kde_n_R, $kde_h_opt) = $cluster1d_obj->one_d_2cluster();
printf("# clustering of   d: %5d  k-means: %5d below %5d above %8.6f, q: %6.4f;  kde: %5d below %5d above %8.6f.\n", $n_pts, $km_n_L, $km_n_R, $km_h_opt, $q, $kde_n_L, $kde_n_R, $kde_h_opt);
my $max_ok_d = $km_h_opt;

# *****  Look at alternative pedigrees if requested
if ($find_alternatives > 0) {
  $max_self_agmr = sprintf("%.6f", $max_self_agmr);
  $max_ok_hgmr = sprintf("%.6f", $max_ok_hgmr);
  $max_self_r = sprintf("%.6f", $max_self_r);
  $max_ok_z = sprintf("%.6f", $max_ok_z);
  $max_ok_d = sprintf("%.6f", $max_ok_d);
    
#  $gtfilename = 'genotype_matrix_out';
  $command = "pedigree_test -g $gtfilename -p $pedigree_table_filename   -x $max_bad_gt_fraction -A $find_alternatives -o $c_output_filename_alt ";
  $command .= " -a $max_self_agmr  -h $max_ok_hgmr -r $max_self_r -D $max_ok_d";
  print "# command: $command \n";
  system "$command";


  open my $fh_alt, "<", "$c_output_filename_alt";
  @lines = <$fh_alt>; # read the pedigrees_and_alternatives file into array

  my $factor = 0.33; # denominators will be considered 'good' if >= factor*median denominator for the quantity
  my %category_counts = ();
  my $bad_denoms_count = 0;
  my %filename_handle = ();

  open my $fhbest, ">", "$best_pedigrees_filename";
  while (my ($j, $line) = each @lines) {
    next if($line =~ /^\s*#/);
    my @cols = split(" ", $line);
    my ($acc_id, $acc_md) = @cols[0,1];
    my ($mat_id, $pat_id) = @cols[2,3];
    my $id_pair = "$mat_id $pat_id";
    my $ped_d = $cols[17];

    my %allped_d = ();
    my $ok_pedigrees_count = 0; # counts all ok (small d) pedigrees for this accession, both pedigree from table and alternatives.
#    my $denoms_ok = are_denoms_ok(\@cols, 4, $factor, $median_matpat_agmr_denom, $median_hgmr_denom, $median_r_denom, $median_z_denom, $median_d_denom);
    my $category_string = category(\@cols, 4, $max_self_agmr, $max_ok_hgmr, $max_self_r, $max_ok_z, $max_ok_d);
    my $no_space_cat_str = $category_string;
    $no_space_cat_str =~ s/ /_/g;
    if(exists $filename_handle{$no_space_cat_str}){
      my $fh = $filename_handle{$no_space_cat_str};
       print $fh join("  ", @cols[0..17]), "\n";
    }else{
      open my $fh, ">", "$no_space_cat_str";
      $filename_handle{$no_space_cat_str} = $fh;
        print $fh join("  ", @cols[0..17]), "\n";
    }
   
    my $ped_str = sprintf("  ped  %20s %20s  ", $mat_id, $pat_id) . "  $category_string";
    if (1){ #  or  $denoms_ok) {
      $ok_pedigrees_count++ if ($category_string eq '0 00 00 00'  or  $category_string eq '1 01 01 00');
      $category_counts{$category_string}++;
    } else {
      $bad_denoms_count++;
   #   print STDERR join("  ", @cols[0..17]), "\n";
      $category_counts{'x xx xx xx'}++;
    }
    $allped_d{$ped_str} = $ped_d if(looks_like_number($ped_d));

    my $n_cols_per_ped = 16; # Fpar id, Mpar id, denom and numer/denom for agmrFM, hgmrF, rF, hgmrM, rM, z, d
    # alternative pedigrees:
    my $n_alternatives = $cols[$n_cols_per_ped+2] // 0; #
    my %okalt_d = ();
    for (my $i = 0; $i < $n_alternatives; $i++) {
      my $first = $n_cols_per_ped+5 + $i*$n_cols_per_ped;
      my $alt_category_string = '';
      my $alt_id_pair = sprintf("%20s %20s", $cols[$first-2],  $cols[$first-1]);
      my $alt_d = $cols[$first+13];
 #     $denoms_ok = are_denoms_ok(\@cols, $first, $factor, $median_matpat_agmr_denom, $median_hgmr_denom, $median_r_denom, $median_z_denom, $median_d_denom);
      if (1){ #  or  $denoms_ok) {
	$alt_category_string = category(\@cols, $first, $max_self_agmr, $max_ok_hgmr, $max_self_r, $max_ok_z, $max_ok_d);
	if ($alt_category_string eq '0 00 00 00'  or $alt_category_string eq '1 01 01 00') {
	  $ok_pedigrees_count++;
	  $okalt_d{"  alt  $alt_id_pair  $alt_category_string"} = $alt_d if(looks_like_number($alt_d));
	}
      } else {
	$alt_category_string = 'x xx xx xx';
      }
      $allped_d{"  alt  $alt_id_pair  $alt_category_string"} = $alt_d if(looks_like_number($alt_d));
    }

    my $output_string = $cols[0] . "  $ok_pedigrees_count  ";
    my @sorted_alts = #sort {$okalt_d{$a} <=> $okalt_d{$b} } keys %okalt_d;
      sort { # sort first by number of x's (low to high), then by d (low_to_high)
	my $acat = ($a =~ /(\S+\s+\S+\s+\S+\s+\S+)\s*$/)? $1 : 'x xx xx xx';
	my $bcat = ($b =~ /(\S+\s+\S+\s+\S+\s+\S+)\s*$/)? $1 : 'x xx xx xx';
	my $ax = () = $acat =~ /x/g;
	my $bx = () = $bcat =~ /x/g;
	(($ax <=> $bx) or ($okalt_d{$a} <=> $okalt_d{$b})); } keys %okalt_d;
    
    my @sorted_allpeds = # sort {$allped_d{$a} <=> $allped_d{$b} } keys %allped_d;
      sort {
	my $acat = ($a =~ /(\S+\s+\S+\s+\S+\s+\S+)\s*$/)? $1 : 'x xx xx xx';
	my $bcat = ($b =~ /(\S+\s+\S+\s+\S+\s+\S+)\s*$/)? $1 : 'x xx xx xx';
	my $ax = () = $acat =~ /x/g;
	my $bx = () = $bcat =~ /x/g;
	(($ax <=> $bx) or ($allped_d{$a} <=> $allped_d{$b}));
      } keys %allped_d;
      
    if ($category_string eq '0 00 00 00'  or  $category_string eq '1 01 01 00') { # pedigree from table is 'good'
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
  my $total_count = 0;
  for my $cat (@scategories) {
    my $cat_count = $category_counts{$cat};
    print "$cat  $cat_count \n";
    $total_count += $cat_count;
  }
  print "total:  $total_count\n";
}


sub are_denoms_ok{
  my @cols = @{my $cls = shift};
  my $first = shift; # the col with the matpat_agmr denom.
  my $factor = shift; # denom is ok if it is >= factor*median_denom
  my ($median_matpat_agmr_denom, $median_hgmr_denom, $median_r_denom, $median_z_denom, $median_d_denom) = @_;
  my $last = $first + 13; # col with d
  my ($FMagmr_denom, $FMagmr, $Fhgmr_denom, $Fhgmr, $Fr_denom, $Fr, $Mhgmr_denom, $Mhgmr, $Mr_denom, $Mr, $z_denom, $z, $d_denom, $d) = @cols[$first..$last];
  my $agmr_denom_ok = ($FMagmr_denom >= $factor*$median_matpat_agmr_denom);
  my $hgmr_denoms_ok = ($Fhgmr_denom >= $factor*$median_hgmr_denom  and  $Mhgmr_denom >= $factor*$median_hgmr_denom);
  my $r_denoms_ok = ($Fr_denom >= $factor*$median_r_denom  and  $Mr_denom >= $factor*$median_r_denom);
  my $z_denom_ok = ($z_denom >= $factor*$median_z_denom);
  my $d_denom_ok = ($d_denom >= $factor*$median_d_denom);
  my $denoms_ok = ($agmr_denom_ok  and  $hgmr_denoms_ok  and  $r_denoms_ok  and  $z_denom_ok  and  $d_denom_ok);
  return $denoms_ok;
}

sub category{
  my @cols = @{my $cls = shift};
  my $first = shift;
  my ($max_self_agmr, $max_ok_hgmr, $max_self_r, $max_ok_z, $max_ok_d) = @_;
  my $last = $first + 13;
  my ($FMagmr_denom, $FMagmr, $Fhgmr_denom, $Fhgmr, $Fr_denom, $Fr, $Mhgmr_denom, $Mhgmr, $Mr_denom, $Mr, $z_denom, $z, $d_denom, $d) = @cols[$first..$last];
  my $category_string = '';
  # $category_string .= ($FMagmr <= $max_self_agmr)? '0' : '1';
  # $category_string .= ' ';
  # $category_string .= ($Fhgmr <= $max_ok_hgmr)? '0' : '1';
  # $category_string .= ($Fr <= $max_self_r)? '0' : '1';
  # $category_string .= ' ';
  # $category_string .= ($Mhgmr <= $max_ok_hgmr)? '0' : '1';
  # $category_string .= ($Mr <= $max_self_r)? '0' : '1';
  # $category_string .= ' ';
  # $category_string .= ($z <= $max_ok_z)? '0' : '1';
  # $category_string .= ($d <= $max_ok_d)? '0' : '1';

  $category_string .= category_character($FMagmr_denom, $FMagmr, $max_self_agmr);
  $category_string .= ' ';
   $category_string .= category_character($Fhgmr_denom, $Fhgmr, $max_ok_hgmr);
   $category_string .= category_character($Fr_denom, $Fr, $max_self_r);
  $category_string .= ' ';
   $category_string .= category_character($Mhgmr_denom, $Mhgmr, $max_ok_hgmr);
   $category_string .= category_character($Mr_denom, $Mr, $max_self_r);
  $category_string .= ' ';
   $category_string .= category_character($z_denom, $z, $max_ok_z);
  $category_string .= category_character($d_denom, $d, $max_ok_d);
  
  return $category_string;
}

sub category_character{
  my $denom = shift;
  my $r = shift;
  my $max = shift;
  if($denom <= 0){
    return 'x';
  }else{
    return ($r <= $max)? '0' : '1';
  }
}

