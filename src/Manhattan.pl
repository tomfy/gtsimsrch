#!/usr/bin/perl -w
use strict;

# make a manhattan plot from a gwas output file with this format:
# 1       S1_84637        84637   H       L       0.0608919       -0.00179986     0.0912604       0.984265
# col 1 is chromosome number, and col 3 is position within the chromosome.
# col 9 is the p value.

my $causal_snp_file = shift // undef;
my $out_file = 'posp';
my $persist = 1;
my $width = 1000;
my $height = 700;

my %causalsnp = ();
if(defined $causal_snp_file){
open my $fh_causal, "<", "$causal_snp_file";
while(<$fh_causal>){
  next if(/^\s*SNP/);
  my @cols = split(" ", $_);
  $causalsnp{$cols[0]} = 1;
}
}

my $n_snps = 0;
my %chrom__pos_p = ();
while(<>){
  next if(/^\s*(Chr|#)/);
  my @cols = split(" ", $_);
  my ($chrom, $id, $pos, $p) = @cols[0,1,2,8];
  if(!exists $chrom__pos_p{$chrom}){
    $chrom__pos_p{$chrom} = {};
  }
  $chrom__pos_p{$chrom}->{$pos} = [$id, $p];
  $n_snps++;
}

my $max_cume_position = 0;
open my $fh_out, ">", "$out_file";
my $chrom_start_position = 0;
my @chroms = sort {$a <=> $b} keys %chrom__pos_p;
for my $achrom (@chroms){
  my %pos_p = %{$chrom__pos_p{$achrom}};
  my @positions = sort {$a <=> $b} keys %pos_p;
  my $max_position_on_chrom = $positions[-1];
  for my $a_pos (@positions){
    my ($id, $p_value) = @{$pos_p{$a_pos}};
    my $causal = (exists $causalsnp{$id})? 1 : 0;
    print $fh_out "$id  $achrom  $a_pos  ", $chrom_start_position + $a_pos, "  $p_value  $causal\n";
  }
  $chrom_start_position += $max_position_on_chrom;
}

system "gnuplot -c ~/gtsimsrch/src/manhattan.gnuplot $width $height $out_file $n_snps $chrom_start_position $persist";

