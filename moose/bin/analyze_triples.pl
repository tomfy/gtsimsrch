#!/usr/bin/perl -w
use strict;
use List::Util qw(min max sum);

# store the info for each parent-progeny pair and each chromosome
# from phased_parents.pl output.

my $pedigree_test_output = shift;

my $max_recomb_rate = shift // 0.18;

my %ProgFMpars_info = ();
#my %Mprog_parent_pairs = ();

open my $fhpt, "<", "$pedigree_test_output";
while (my $pt_line = <$fhpt>) {
  chomp $pt_line;
  my @cols = split(" ", $pt_line);
  my ($A, $F, $M) = @cols[0, 2, 3];
   $ProgFMpars_info{"$A $F $M"} .= "  $pt_line";
  # $Mprog_parent_pairs{"$A $M"} .= "  $pt_line";
}
close $fhpt;


my %pppair_info = ();
while (<>) {
  my @cols = split(" ", $_);
  my ($prog, $par, $chrom, $XA, $XB, $Nhet) = @cols[2,3,4,6,7,10];
  my $ppp = "$par $prog";
  # next if(min($XA, $XB)/$Nhet > $max_recomb_rate); # donot store if poor parent-progeny candidate
  my $the_info = {XA => $XA, XB => $XB, Nhet => $Nhet};
  # print "$XA $XB  $Nhet\n";
  if (!exists $pppair_info{$ppp}) {
    $pppair_info{$ppp} = {};
  }
  $pppair_info{$ppp}->{$chrom} = $the_info;
}
print STDERR "done with input of ", scalar keys %pppair_info, " parent-progeny pairs.\n";


my @pppairs = keys %pppair_info;
for my $ppp (@pppairs) {
  my $chr_info = $pppair_info{$ppp};
  my ($Xmin_sum, $Nhet_sum) = (0, 0);
  while (my($chr, $info) = each %$chr_info) {
    $Xmin_sum += min($info->{XA}, $info->{XB});
    $Nhet_sum += $info->{Nhet};
  }
#  print STDERR "$ppp   $Xmin_sum $Nhet_sum  ", $Xmin_sum/$Nhet_sum, "\n";
  delete($pppair_info{$ppp}) if(($Nhet_sum == 0)  or  $Xmin_sum/$Nhet_sum  > $max_recomb_rate);
}
#exit;
@pppairs = keys %pppair_info;
print STDERR "### N pp pairs kept: ", scalar @pppairs, "\n";
while (my($ippp1, $pppair1) = each @pppairs) {
  if ($pppair1 =~ /^\s*(\S+)\s+(\S+)\s*$/) {
    my ($parent1, $progeny) = ($1, $2);
  #  for (my $ippp2 = $ippp1; $ippp2 < scalar @pppairs; $ippp2++) {
       for (my $ippp2 = 0; $ippp2 < scalar @pppairs; $ippp2++) {
      my $pppair2 = $pppairs[$ippp2];
      if ($pppair2 =~ /^\s*(\S+)\s+(\S+)\s*$/) {
	my ($parent2, $progeny2) = ($1, $2);
	next if($progeny2 ne $progeny);
	my $key = "$progeny $parent1 $parent2";
	my $PFMi = $ProgFMpars_info{$key} // 'XXX';
	next if($PFMi eq 'XXX');
	# get figure of merit for this triple:
	my ($Xmin1, $Xmin2, $Xbest, $Nhet1, $Nhet2, $Nbad_chrom) = triple_quality(\%pppair_info, $parent1, $parent2, $progeny);
	my $Nhet_total = $Nhet1 + $Nhet2;
	print "$progeny   $parent1 $parent2   $Xmin1 $Xmin2 $Xbest  $Nhet1 $Nhet2  ", $Xmin1/$Nhet1, "  ", $Xmin2/$Nhet2, "  ", $Xbest/$Nhet_total, "  $Nbad_chrom  $PFMi \n";
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
   my $Xmin1 = 0;
  my $Xmin2 = 0;
  my $Xbest = 0;
  my $Xmin1_sum = 0;
  my $Xmin2_sum = 0;
  my $Xbest_sum = 0;
  my $Nhet1 = 0;
  my $Nhet2 = 0;
  my $bad_count = 0;
  for my $i_chrom (@chroms){
    # print STDERR "$i_chrom  ", join(" ", keys %{$chrom_info1->{$i_chrom}}), "  ",  $chrom_info1->{$i_chrom}->{XA}, "  ", $chrom_info2->{$i_chrom}->{XB}, "\n";
    $Xmin1 = min( $chrom_info1->{$i_chrom}->{XA}, $chrom_info1->{$i_chrom}->{XB} );
    $Xmin2 = min( $chrom_info2->{$i_chrom}->{XA}, $chrom_info2->{$i_chrom}->{XB} );
    $Xmin1_sum += $Xmin1;
    $Xmin2_sum += $Xmin2;
    my $X_AB = $chrom_info1->{$i_chrom}->{XA} + $chrom_info2->{$i_chrom}->{XB};
    my $X_BA = $chrom_info1->{$i_chrom}->{XB} + $chrom_info2->{$i_chrom}->{XA};
    $Xbest = min($X_AB, $X_BA);
    $Xbest_sum += $Xbest;
    $Nhet1 += $chrom_info1->{$i_chrom}->{Nhet};
    $Nhet2 += $chrom_info2->{$i_chrom}->{Nhet};
    $bad_count++ if($Xmin1+$Xmin2 < $Xbest);
  }
  return ($Xmin1_sum, $Xmin2_sum, $Xbest_sum, $Nhet1, $Nhet2, $bad_count);
}
