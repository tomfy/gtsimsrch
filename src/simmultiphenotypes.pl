#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw(min max sum shuffle);

my $binary_fileset = undef;
my $n_phenotypes = 2;
my $causal_snps_file = undef; # if not specified, choose $n_causal_snps at random to be causal.
my $pheno_out_filename = undef;
my $verbose = 0;
my $effect_size = 0.1;
GetOptions(
	   'bin_fileset|fileset_bin=s' => \$binary_fileset,
	   'causal_snps_file=s' => \$causal_snps_file,
	   'n_phenotypes=i' => \$n_phenotypes,
	   'pheno_out_filename=s' => \$pheno_out_filename,
	   'verbose!' => \$verbose,
	   'effect_size=s' => \$effect_size,
	  );

my %accid_outstring = ();
my $phenotype_count = 0;
# from the full genotype set (with duplicates) simulate data with $n_causal_snps
my %accid_phenotypes = ();
my $outfile = $binary_fileset;

my $causal_variants_filename;
if (!defined $causal_snps_file) {
  die;
} else {
  $causal_variants_filename = $causal_snps_file . "_id_effect";
  fix_causal_variants_file($causal_snps_file, $effect_size, $causal_variants_filename);
}

# print STDERR `cat $causal_variants_filename `, "\n"; # getc();

my $ref_accid_phenotype = {};
my $accid_phenotype = {};
for my $i (1..$n_phenotypes) {
  my $phenotype_filename = # $i . "_" .
    $outfile . ".pheno";
  my $outfile_i = # $i . "_" .
    $outfile;
  my $simu_command = "Simu  --bfile $binary_fileset  --qt  --causal-variants $causal_variants_filename  --out $outfile_i ";
  my $simu_stdout = `$simu_command`;
  if ($i == 1) {
    $ref_accid_phenotype = phenotypes_from_file($phenotype_filename);
    #print STDERR "size 1: ", scalar keys %$ref_accid_phenotype, "\n";
    add_phenotypes($ref_accid_phenotype, \%accid_phenotypes);
  } else {
    $accid_phenotype = phenotypes_from_file($phenotype_filename);
    #print STDERR "size 2: ", scalar keys %$accid_phenotype, "\n";
    my $the_slope = slope($ref_accid_phenotype, $accid_phenotype);
    print STDERR "slope: before: $the_slope;   ";
    if ($the_slope < 0) {
      flip_signs($accid_phenotype);
    }
    add_phenotypes($accid_phenotype, \%accid_phenotypes);
    $the_slope = slope($ref_accid_phenotype, $accid_phenotype);
    print STDERR "after: $the_slope\n";
  }
  $phenotype_count++;
}
die if($phenotype_count !=  $n_phenotypes);
while (my ($id, $phs) = each %accid_phenotypes) {
  die if(scalar @$phs !=  $n_phenotypes);
  # print STDERR "$id  $id " . sum(@$phs)/$n_phenotypes . "\n";
  $accid_outstring{$id} = "$id  $id " . sum(@$phs)/$n_phenotypes;
  if ($verbose) {
    $accid_outstring{$id} .= "   $n_phenotypes  " . join(" ", @$phs);
  }
}

my @sorted_ids = sort keys %accid_outstring;
#print STDERR "will output ", scalar @sorted_ids, "  phenotypes.\n";
$pheno_out_filename = $n_phenotypes . "_phenotypes" if(!defined $pheno_out_filename);
open my $fhout, ">", "$pheno_out_filename" or die "Couldn't open $pheno_out_filename for writing. \n";
for my $anid (@sorted_ids) {
  print $fhout  $accid_outstring{$anid}, "\n";
}
close $fhout;


sub add_phenotypes{
  my $accid_ph = shift;
  my $accid_phenos = shift;
  while (my ($id, $pheno) = each %$accid_ph) {
    $accid_phenos->{$id} = [] if(!exists $accid_phenos->{$id});
    push @{ $accid_phenos->{$id} }, $pheno;
  }
}

sub fix_causal_variants_file{
  my $infile = shift;
  my $effect_size = shift;
  my $outfile = shift;
  open my $fh, "<", "$infile";
  open my $fhout, ">", "$outfile";
  while (<$fh>) {
    next if(/^SNP/);
    my @cols = split(" ", $_);
    print STDERR "effect size: $effect_size \n";
    if($effect_size ne 'R'){
      print $fhout $cols[0], "  $effect_size\n";
      $effect_size *= -0.5;
    }else{
      print $fhout $cols[0], "\n";
    }
  }
  close $fh;
  close $fhout;
}

sub phenotypes_from_file{
  my $pheno_file = shift;
  my %accid_phenotype = ();
  open my $fh, "<", "$pheno_file";
  while (<$fh>) {
    next if(/^\s*FID/);
    my ($fid, $iid, $pheno) = split(" ", $_);
    $accid_phenotype{$fid} = $pheno;
    # print STDERR "$fid $pheno\n";
  }
  close $fh;
  return \%accid_phenotype;
}

sub slope{
  my $aph1 = shift;
  my $aph2 = shift;
  my $size = scalar keys %$aph1;
  my $size2 = scalar keys %$aph2;
  die "$size $size2\n" if($size != $size2);
  my $avg_phe1 = sum(values %$aph1)/$size;
  my $avg_phe2 = sum(values %$aph2)/$size;
  my $numerator = 0;
  my $denominator = 0;
  while (my($id, $aphe1) = each %$aph1) {
    my $aphe2 = $aph2->{$id};
    $numerator += ($aphe1 - $avg_phe1)*($aphe2 - $avg_phe2);
    $denominator += ($aphe1 - $avg_phe1)**2;
  }
  return $numerator/$denominator;
}

sub flip_signs{ # given a array ref (of numbers), flip the sign of each element.
  my $hashref = shift;
  while (my($k, $v) = each %$hashref) {
    $hashref->{$k} *= -1;
  }
}



