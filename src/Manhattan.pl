#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw(min max sum shuffle);
use File::Spec qw(splitpath);
use File::Basename 'dirname';

use Cwd 'abs_path';
my $bindir;
BEGIN {     # this has to go in Begin block so happens at compile time
  $bindir =
    dirname( abs_path(__FILE__) ) ; # the directory containing this script
}

# make a manhattan plot from a gwas output file with this format:
# 1       S1_84637        84637   H       L       0.0608919       -0.00179986     0.0912604       0.984265
# col 1 is chromosome number, and col 3 is position within the chromosome.
# col 9 is the p value.
# if   -causal_snps_file  <snp id filename>
# will make the snps listed in snp_id_filename a different color.

my $causal_snps_file = undef;
my $temp_out_file = 'temp';
my $persist = 1;
my $width = 800;
my $height = 500;
my $terminal = 'x11';
my $png_filename = "mnhttn.png";
my $title = undef;
my $gwas_output_file = undef;

GetOptions(
	   'gwas_out_file=s' => \$gwas_output_file,
	   'causal_snps_file=s' => \$causal_snps_file,

	   'png_filename=s' => \$png_filename,
	   'title=s' => \$title,

	   'width=i' => \$width,
	   'height=i' => \$height,
	   'terminal=s' => \$terminal,

	   'temp_output_file=s' => \$temp_out_file,
	   );

my %causalsnp = ();
if(defined $causal_snps_file){
open my $fh_causal, "<", "$causal_snps_file";
while(<$fh_causal>){
  next if(/^\s*SNP/);
  my @cols = split(" ", $_);
  $causalsnp{$cols[0]} = 1;
}
}

my $n_snps = 0;
my %chrom__pos_p = ();
open my $fh_gwas, "<", "$gwas_output_file";
while(<$fh_gwas>){
  next if(/^\s*(C[hH][rR]|#)/);
  my @cols = split(" ", $_);
  my ($chrom, $id, $pos, $p) = @cols[0,1,2,8];
  if(!exists $chrom__pos_p{$chrom}){
    $chrom__pos_p{$chrom} = {};
  }
  $chrom__pos_p{$chrom}->{$pos} = [$id, $p];
  $n_snps++;
}

my $max_cume_position = 0;
open my $fh_out, ">", "$temp_out_file";
my $chrom_start_position = 0;
my @chroms = sort {$a <=> $b} keys %chrom__pos_p;
for my $achrom (@chroms){
  my %pos_p = %{$chrom__pos_p{$achrom}};
  my @positions = sort {$a <=> $b} keys %pos_p;
  my $max_position_on_chrom = $positions[-1];
  for my $a_pos (@positions){
    my ($id, $p_value) = @{$pos_p{$a_pos}};
    my $causal = 0;
    $causal = 1 if(defined $causal_snps_file  and  exists $causalsnp{$id});
    print $fh_out "$id  $achrom  $a_pos  ", $chrom_start_position + $a_pos, "  $p_value  $causal\n";
  }
  $chrom_start_position += $max_position_on_chrom;
}

if(!defined $title){
  $title = $gwas_output_file;
}

my $manhattan_path = $bindir . '/manhattan.gnuplot ';
system "gnuplot -c $manhattan_path  $width $height $temp_out_file $n_snps $chrom_start_position $persist $terminal $png_filename $title";
