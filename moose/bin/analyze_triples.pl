#!/usr/bin/perl -w
use strict;
use List::Util qw(min max sum);

# store the info for each parent-progeny pair and each chromosome
# from phased_parents.pl output.

my $max_recomb_rate = shift // 0.125;

my %pppair_info = ();
while (<>) {
  my @cols = split(" ", $_);
  my ($par, $prog, $chrom, $XA, $XB, $Nhet) = @cols[2,3,4,6,7,10];
  my $ppp = "$par $prog";
  # next if(min($XA, $XB)/$Nhet > $max_recomb_rate); # do not store if poor parent-progeny candidate
  my $the_info = {XA => $XA, XB => $XB, Nhet => $Nhet};
  # print "$XA $XB  $Nhet\n";
  if (!exists $pppair_info{$ppp}) {
    $pppair_info{$ppp} = {};
  }
  $pppair_info{$ppp}->{$chrom} = $the_info;
}



my @pppairs = keys %pppair_info;
for my $ppp (@pppairs) {
  my $chr_info = $pppair_info{$ppp};
  my ($Xmin_sum, $Nhet_sum) = (0, 0);
  while (my($chr, $info) = each %$chr_info) {
    $Xmin_sum += min($info->{XA}, $info->{XB});
    $Nhet_sum += $info->{Nhet};
  }
#  print STDERR "$ppp   $Xmin_sum $Nhet_sum  ", $Xmin_sum/$Nhet_sum, "\n";
  delete($pppair_info{$ppp}) if($Xmin_sum/$Nhet_sum  > $max_recomb_rate);
}
#exit;
@pppairs = keys %pppair_info;

while (my($ippp1, $pppair1) = each @pppairs) {
  if ($pppair1 =~ /^\s*(\S+)\s+(\S+)\s*$/) {
    my ($parent1, $progeny) = ($1, $2);
    for (my $ippp2 = $ippp1; $ippp2 < scalar @pppairs; $ippp2++) {
      my $pppair2 = $pppairs[$ippp2];
      if ($pppair2 =~ /^\s*(\S+)\s+(\S+)\s*$/) {
	my ($parent2, $progeny2) = ($1, $2);
	next if($progeny2 ne $progeny);
	# get figure of merit for this triple:
	my ($Xbest, $Nhet_total) = triple_quality(\%pppair_info, $parent1, $parent2, $progeny);
	print "$progeny   $parent1 $parent2   $Xbest  $Nhet_total  ", $Xbest/$Nhet_total, "\n";
      }
    }
  }
}



sub triple_quality{
  my $ppp_info = shift;
  my $par1 = shift;
  my $par2 = shift;
  my $prog = shift;
  my $chrom_info1 = $ppp_info->{"$par1 $prog"};
  my $chrom_info2 = $ppp_info->{"$par2 $prog"};
  my @chroms = sort {$a <=> $b} keys %$chrom_info1;
  my $Xbest = 0;
  my $Nhet_total = 0;
  for my $i_chrom (@chroms){
    print STDERR "$i_chrom  ", join(" ", keys %{$chrom_info1->{$i_chrom}}), "  ",  $chrom_info1->{$i_chrom}->{XA}, "  ", $chrom_info2->{$i_chrom}->{XB}, "\n";
    my $X_AB = $chrom_info1->{$i_chrom}->{XA} + $chrom_info2->{$i_chrom}->{XB};
    my $X_BA = $chrom_info1->{$i_chrom}->{XB} + $chrom_info2->{$i_chrom}->{XA};
    $Xbest += min($X_AB, $X_BA);
    $Nhet_total += $chrom_info1->{$i_chrom}->{Nhet} + $chrom_info2->{$i_chrom}->{Nhet};
  }
  return ($Xbest, $Nhet_total);
}
