#!/usr/bin/perl -w
use strict;
use List::Util qw(min max sum);
use Getopt::Long;

# store the info for each parent-progeny pair and each chromosome
# from phased_parents.pl output.

my $pedigree_test_output = undef;
my $pp_file = undef;
my $output_file = "antrip.out";
my $max_crossover_rate = 100;
my $help = 0;


GetOptions(
	   'pedigree_test|pt=s' => \$pedigree_test_output,
	   'pp_out|phased_parents_out=s' => \$pp_file,
	   'output_file=s' => \$output_file,
	   'max_crossover_rate=f' => \$max_crossover_rate,
	   'help!' => \$help,
	   );

if(!defined $pp_file  or  !defined $pedigree_test_output  or  $help){
  print STDOUT "usage: analyze_triples.pl -ped <pedigree_test output file>  -pp_out <phased_parents output file> [-out <output file>] \n";
  exit();
}

# ###  read and store pedigree_test output #############
my %ProgFMpars_info = ();
open my $fhpt, "<", "$pedigree_test_output";
while (my $pt_line = <$fhpt>) {
  chomp $pt_line;
  my @cols = split(" ", $pt_line);
  my ($A, $F, $M) = @cols[0, 2, 3];
  # print "$A $F $M\n";
  $ProgFMpars_info{"$A $F $M"} .= "  $pt_line";
  #print STDERR "$A  $F $M\n";
}
close $fhpt;
print STDERR "Pedigrees stored: ", scalar keys %ProgFMpars_info, "\n";
# #####################################################

# ###  read and store output from phased_parents.pl ################
open my $fhpp, "<", "$pp_file";
my %pppair_info = ();
while (my $pp_line = <$fhpp>) {
  next if($pp_line =~ /^\s*#/);
  my @cols = split(" ", $pp_line);
  my ($prog, $par, $chrom, $XA, $XB, $Nhet, $type, $for_rev) = @cols[0,1,2,4,5,8,9,10]; # @cols[2,3,4,6,7,10];
  my $ppp = "$par $prog";
  # next if(min($XA, $XB)/$Nhet > $max_recomb_rate); # donot store if poor parent-progeny candidate
  my $the_info = {XA => $XA, XB => $XB, Nhet => $Nhet};
  # print "$XA $XB  $Nhet\n";
  if (!exists $pppair_info{$ppp}) {
    $pppair_info{$ppp} = {};
  }
  $pppair_info{$ppp}->{$chrom} = $the_info;
}
close $fhpp;
print STDERR "done with input of ", scalar keys %pppair_info, " parent-progeny pairs.\n";
####################################################################

# ###  
my @pppairs = keys %pppair_info; # phased_parents info
for my $ppp (@pppairs) {
  my $chr_info = $pppair_info{$ppp};
  my ($Xmin_sum, $Nhet_sum) = (0, 0);
  #print STDERR "size of chr_info: ", scalar %$chr_info, "\n";
  while (my($chr, $info) = each %$chr_info) {
    $Xmin_sum += min($info->{XA}, $info->{XB});
    $Nhet_sum += $info->{Nhet};
  }
#  print STDERR "$ppp   $Xmin_sum $Nhet_sum\n";
  if (($Nhet_sum == 0)  or  $Xmin_sum/$Nhet_sum  > $max_crossover_rate) {
    print STDERR "ppp:  $ppp   Xmin_sum: $Xmin_sum  Nhet_sum: $Nhet_sum    \n";
    delete($pppair_info{$ppp});
  }else{
 # print STDERR "$ppp   $Xmin_sum $Nhet_sum  ", $Xmin_sum/$Nhet_sum, "\n";
  }
}
#exit;
@pppairs = keys %pppair_info; 
print STDERR "### N pp accession pairs kept: ", scalar @pppairs, "\n";
# exit;

  open my $fhout, ">", "$output_file";
  while (my($AFM, $info) = each %ProgFMpars_info) {
    my ($A, $F, $M) = split(" ", $AFM);
    # my $AF_ppinfo = $pppair_info{"$A $F"};
    #print STDERR "\n\n"; print STDERR "AFM:  $A  $F  $M\n";
    my ($Xmin1, $Xmin2, $Xb1, $Xb2, $Xbest, $Nhet1, $Nhet2, $Nbad_chrom) = triple_quality(\%pppair_info, $F, $M, $A);
    # print STDERR "   $Xmin1 $Xmin2  $Xbest  $Nhet1  $Nhet2  $Nbad_chrom\n";
    my $Nhet_total = $Nhet1 + $Nhet2;
    if ($Nhet1 > 0  and $Nhet2 > 0) {
      print $fhout "$A $F $M  $Xmin1  $Xmin2  " .
	  #  $Xb1 $Xb2
	  "$Xbest  " .
	" $Nhet1 $Nhet2  ", $Xmin1/$Nhet1, "  ", $Xmin2/$Nhet2, "  ", $Xbest/$Nhet_total, "  $Nbad_chrom  $info ", ($F eq $M)? 'S' : 'B', "\n";
    }
  }
close $fhout;



sub triple_quality{
  my $ppp_info = shift;
  my $par1 = shift;
  my $par2 = shift;
  my $prog = shift;
  my $chrom_info1 = $ppp_info->{"$par1 $prog"};
  my $chrom_info2 = $ppp_info->{"$par2 $prog"};
  #print STDERR "par2 prog:   $par2 $prog\n";
 
  my @chroms = sort {$a <=> $b} keys %$chrom_info1;
  my $Xmin1 = 0;
  my $Xmin2 = 0;
  my $Xbest = 0;
  my ($Xb1, $Xb2) = (0, 0); # numbers of Xovers for parent1, parent2 for best solution for triple. ($Xb1 + $Xb2 == $Xbest_sum)
  my $Xmin1_sum = 0;
  my $Xmin2_sum = 0;
  my $Xbest_sum = 0;
  my $Nhet1 = 0;
  my $Nhet2 = 0;
  my $bad_count = 0;
  for my $i_chrom (@chroms) {
    #print STDERR "$i_chrom  ", join(" ", keys %{$chrom_info1->{$i_chrom}}), "  ",  $chrom_info1->{$i_chrom}->{XA} // '-', "  ", $chrom_info2->{$i_chrom}->{XB} // '-', "\n";
    # print STDERR "chrom_info2: ", $chrom_info2
    if(exists $chrom_info1->{$i_chrom}->{XA}  and  exists $chrom_info1->{$i_chrom}->{XB}  and exists $chrom_info2->{$i_chrom}->{XA}  and  exists $chrom_info2->{$i_chrom}->{XB}){
    $Xmin1 = min( $chrom_info1->{$i_chrom}->{XA}, $chrom_info1->{$i_chrom}->{XB} );
    $Xmin2 = min( $chrom_info2->{$i_chrom}->{XA}, $chrom_info2->{$i_chrom}->{XB} );
  }else{
    return (-1, -1, -1, -1, -1, -1);
  }
    $Xmin1_sum += $Xmin1;
    $Xmin2_sum += $Xmin2;
    my $X_AB = $chrom_info1->{$i_chrom}->{XA} + $chrom_info2->{$i_chrom}->{XB};
    my $X_BA = $chrom_info1->{$i_chrom}->{XB} + $chrom_info2->{$i_chrom}->{XA};
    if($X_AB < $X_BA){
      $Xbest = $X_AB;
      $Xb1 += $chrom_info1->{$i_chrom}->{XA};
      $Xb2 += $chrom_info2->{$i_chrom}->{XB}
    }else{
      $Xbest = $X_BA;
      $Xb1 += $chrom_info1->{$i_chrom}->{XB};
      $Xb2 += $chrom_info2->{$i_chrom}->{XA}
    }
    # $Xbest = min($X_AB, $X_BA);
    $Xbest_sum += $Xbest;
    $Nhet1 += $chrom_info1->{$i_chrom}->{Nhet};
    $Nhet2 += $chrom_info2->{$i_chrom}->{Nhet};
    $bad_count++ if($Xmin1+$Xmin2 < $Xbest);
  }
  die "Xbest_sum: $Xbest_sum  Xb1: $Xb1  Xb2: $Xb2\n" if($Xbest_sum != $Xb1 + $Xb2);
  return ($Xmin1_sum, $Xmin2_sum, $Xb1, $Xb2, $Xbest_sum, $Nhet1, $Nhet2, $bad_count);
}
