#!/usr/bin/perl -w
use strict;
use List::Util qw(min max sum shuffle);

# read in a phased dosages file
# i.e. dosage are 0, and 2 for homozygs, and +1, -1 for 2 possible phases of heterozygs.
# analyze to find likely parent-progeny pairs.

use warnings;
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
print STDOUT "# libdir: $libdir \n";
use Chromosome;


my $pedigree_file = undef;
my $pdsgs_file = undef;
my $output_file = "pp.out";
my $rand_parents_to_do = 0; # for each accession with pedigree, also choose this many other pairs of accessions at random to test as parents.
my $do_reverse = 0;
my $max_markers = 5000000000; # to limit to a smaller number of markers for speed in testing
my $use_pedigrees = 0;
my $help = 0;
my $remove_length_one = 0;

GetOptions(
	   'pedigree_file=s' => \$pedigree_file,
	   'pdsgs_file=s' => \$pdsgs_file,
	   'output_file=s' => \$output_file,
	   'rand_parents|n_rand_parents=i' => \$rand_parents_to_do,
	#   'min_prob|gt_min_prob=f' => \$min_gt_prob,
	   'reverse!' => \$do_reverse,
	   'marker_limit=i' => \$max_markers,
	   'help!' => \$help,
	   'remove_length_one!' => \$remove_length_one,
	  );

if (!defined $pedigree_file  or  !defined $pdsgs_file  or  $help) {
  print_help();
  exit;
}

open my $fhout, ">", "$output_file";
my $run_parameter_info = "# pedigree file: $pedigree_file\n".
  "# phased dosages file: $pdsgs_file\n".
  "# output file: $output_file\n".
  "# rand parents per pedigree: $rand_parents_to_do\n".
#  "# min gt probability: $min_gt_prob\n".
  "# do reverse? $do_reverse\n";
print STDOUT $run_parameter_info; print $fhout $run_parameter_info;

##########################################################
# read pedigree file and store accession-parent pairs ####
##########################################################
my %A_Spar = ();
my %A_Fpar = ();
my %A_Mpar = ();
if (defined $pedigree_file) {
  open my $fhped, "<", "$pedigree_file";
  while (<$fhped>) {
    my @cols = split(" ", $_);
    my ($A, $Fpar, $Mpar) = @cols[0,1,2]; # [-3,-2,-1];
    #print STDERR "$A  $Fpar $Mpar\n";
    if ($Fpar ne 'NA'  or  $Mpar ne 'NA') {
      if ($Fpar eq $Mpar) {	# self case
	$A_Spar{$A} = $Fpar;
      } else {			# distinct parents
	$A_Fpar{$A} = $Fpar if($Fpar ne 'NA');
	$A_Mpar{$A} = $Mpar if($Mpar ne 'NA');
      }
    }
    $use_pedigrees = 1;
  }
  close $fhped;
}
print STDOUT  "# Pedigree information from file  $pedigree_file  stored.\n";
print STDOUT "# ", scalar keys %A_Fpar, " female parents;  ", scalar keys %A_Mpar, " male parents.\n";
#sleep(2);
##########################################################
# Done reading pedigree file ######################## ####
##########################################################

##############################################
###  read and store dsgs file  ##############
##############################################

my @acc_ids = ();		# acc ids in pdsgs file
my %acc_chroms = (); # keys: acc ids, values: arrayref of chromosome object;
my %chrom_numbers = ();
#my %acc_gts = ();		# key: acc ids; value:
my @markers = ();
my @chromosome_numbers = (); 

open my $fhpdsgs, "<", "$pdsgs_file";
while (my $line = <$fhpdsgs>) { # get initial 2 lines 
  next if($line =~ /^\s*#/);
  if ($line =~ /^\s*MARKER/) {
    @markers = split(" ", $line);
    shift @markers;
  } elsif ($line =~ /^\s*CHROMOSOME/) {
    @chromosome_numbers = split(" ", $line);
    shift @chromosome_numbers;
    last; 
  }
}
print STDERR "n markers: ", scalar @markers, "  ", scalar @chromosome_numbers, "\n";
for my $ichr (@chromosome_numbers){
  $chrom_numbers{$ichr}++;
}
my $accs_read = 0;
while (my $line = <$fhpdsgs>) # each pass through this loop processes one accessions, all markers.
  {
    my @pdosages = split(" ", $line);
    my $acc_id = shift @pdosages;
    # print STDERR join(" ", @chromosome_numbers), "\n";
    # print STDERR join(" ", @pdosages), "\n";

    # print STDERR scalar @chromosome_numbers, "  ", scalar @pdosages, "\n";
    # print STDOUT "acc_id: $acc_id\n";
    push @acc_ids, $acc_id;
    die if(scalar @chromosome_numbers != scalar @pdosages);
    while (my($i, $pd) = each @pdosages) {
      my $i_chrom = $chromosome_numbers[$i]; # the number of the chromosome which this marker is on.
      if (!exists $acc_chroms{$acc_id}->[$i_chrom]) { # check whether Chromosome obj. exists for this accession and chromosome number
	$acc_chroms{$acc_id}->[$i_chrom] = Chromosome->new({i_chrom => $i_chrom, genotypes => []});
      }
      $acc_chroms{$acc_id}->[$i_chrom]->add_genotype($pd);
    }
    $accs_read++;
    print "Accessions read so far: $accs_read\n" if($accs_read % 100  == 0);
    # for my $ichr (keys %chrom_numbers){
    #   my $chr = $acc_chroms{$acc_id}->[$ichr];
    #   $chr->printxxx();
    # }
    
  }
print STDOUT "# Done storing genotypes for ", scalar @acc_ids, " accessions. \n";
close $fhpdsgs;

###########################################
# done storing info in pdsgs file. ########
###########################################


my @chroms = sort {$a <=> $b} keys %chrom_numbers;
my $n_chrom_pairs_analyzed_forward = 0;

print STDOUT "# use pedigrees? ", ($use_pedigrees)? 'Y' : 'N', "\n";;
while (my ($i, $progid) = each @acc_ids) { # loop over accessions considered as progeny
  my $ped_Fpar;
  my $ped_Mpar;
  if (exists $A_Spar{$progid}) { # self
    $ped_Fpar = $A_Spar{$progid};
    $ped_Mpar = $ped_Fpar;
  } else {
    $ped_Fpar = $A_Fpar{$progid} // 'unknown';
    $ped_Mpar = $A_Mpar{$progid} // 'unknown';
  }
  if ($use_pedigrees  and  ($ped_Fpar eq 'unknown'  and  $ped_Mpar eq 'unknown')) {
    next;
  }
  my %cand_parent_ids = ();
  if ($ped_Fpar eq $ped_Mpar) { # self
    $cand_parent_ids{$ped_Fpar} = 'S';
  } else {
    if ($ped_Fpar ne 'unknown') {
      $cand_parent_ids{$ped_Fpar} = 'F';
    }
    if ($ped_Mpar ne 'unknown') {
      $cand_parent_ids{$ped_Mpar} = 'M';
    }
  }
  my $rand_parents_added = 0;
  my @shuffled_ids = shuffle(@acc_ids);
  for my $anid (@shuffled_ids) { # add some other randomly chosen parents
    last if($rand_parents_added >= $rand_parents_to_do);
    if (exists $acc_chroms{$anid}) {
      if ($anid ne $progid  and  $anid ne $ped_Fpar  and  $anid ne $ped_Mpar) {
	$cand_parent_ids{$anid} = 'R';
	$rand_parents_added++;
      }
    }
  }
  for my $i_chrom (@chroms) {	# loop over chromosomes
    my $pgts1 = $acc_chroms{$progid}->[$i_chrom]->genotypes();
    while (my($parid, $type) = each %cand_parent_ids) {

      next if($parid eq $progid);
      #my $pFM = $type;		# $pF . $pM;
      next if(!exists $acc_chroms{$parid});
      # if ($pFM ne 'XX') {
      #print STDERR "$i_chrom  $progid  $pFM   $ped_Fpar  $ped_Mpar \n";
      # forward direction - $parid is parent
      my $pgts2 = $acc_chroms{$parid}->[$i_chrom]->genotypes();
    
      my ($Do, $So, $XA, $XB, $parent_het_count, $ndubA, $interXlsA, $ndubB, $interXlsB) = analyze_pgts_pair($pgts2, $pgts1); # $pgts1: parent, $pgts2: progeny
      
      #print STDERR "$progid  $parid  $i_chrom  $parent_het_count\n";
      #  print STDERR "  ", join(" ", @$pgts2), "\n";
      # my $L1A = count_L1s($interXlsA);
      # my $L1B = count_L1s($interXlsB);
      #  my $XAprime = $XA - $ndubA;
      #  my $XBprime = $XB - $ndubB;
      if ($remove_length_one) {
	$XA -= $ndubA;
	$XB -= $ndubB;
      }
      $n_chrom_pairs_analyzed_forward++;
      print STDOUT "# Chromosome parent-progeny pairs analyzed: $n_chrom_pairs_analyzed_forward \n" if($n_chrom_pairs_analyzed_forward % 1000 == 0);
      my $hgmr_denom = $Do + $So;
      print $fhout "$progid  $parid  $i_chrom  ",
	 ($hgmr_denom > 0)? $Do/$hgmr_denom : '-1',
	# "$Do $So  ",
	"   $XA $XB   ";
      if ($XA < $XB) {
	print $fhout "$XA $XB   ";
	# while (my($i, $iXl) = each @$interXlsA) {
	#   print STDOUT "A  $iXl  min forward\n";
	# }
	# while (my($i, $iXl) = each @$interXlsB) {
	#   print STDOUT "B  $iXl  max forward\n";
	# }
      } else {			# $XB >= $XA
	print $fhout "$XB $XA   ";
	# while (my($i, $iXl) = each @$interXlsA) {
	#   print STDOUT "A  $iXl  max forward\n";
	# }
	# while (my($i, $iXl) = each @$interXlsB) {
	#   print STDOUT "B  $iXl  min forward\n";
	# }
      }

      #   print $fhout min($XAprime, $XBprime), "  ", max($XAprime, $XBprime), "  ";
      print $fhout " $parent_het_count  $type  forward\n";

      if ($do_reverse) { # now with $progid as parent, $parid as progeny
	my ($Do_rev, $So_rev, $XA_rev, $XB_rev, $parent_het_count_rev,
	    $ndubA,
	    $interXlsA_rev, $ndubB,
	    $interXlsB_rev) = analyze_pgts_pair($pgts1, $pgts2);
	my $hgmr_rev = ($So_rev > 0)? $Do_rev/($Do_rev + $So_rev) : -1;
	if ($remove_length_one) {
	  $XA_rev -= $ndubA;
	  $XB_rev -= $ndubB;
	}
	print $fhout "$parid  $progid  $i_chrom  $hgmr_rev  $XA_rev $XB_rev   ";
	if ($XA_rev < $XB_rev) {
	  print $fhout "$XA_rev $XB_rev  "; # $length1count_A $length1count_B  ";
	} else {
	  print $fhout "$XB_rev $XA_rev  "; #  $length1count_B $length1count_A  ";
	}
	print $fhout " $parent_het_count_rev reverse\n";
      }
    }
  }
}

#  ****************************************************************************

sub analyze_pgts_pair{ # consider $pgts1 as parent, $pgts2 as progeny.
  my $pgts1 = shift;
  my $pgts2 = shift;
  my $chromosome_numbers = shift; 
  die if(scalar @$pgts1 != scalar @$pgts2);
  #  print "par :  ", join('', @$pgts1), "\n";
  #  print "prog:  ", join('', @$pgts2), "\n";
  my $Shom = 0;
  my $Dhom = 0;
  my $phaseA = undef;
  my $phaseB = undef;
  my $previous_phaseA = undef;
  my $previous_phaseB = undef;
  my $switch_countA = 0;
  my $switch_countB = 0;
  my $par_1count = 0;
  my $par_2count = 0;
  my $par_het_count = 0;
  my $length_1_countA = 0;
  my $length_1_countB = 0;
  my ($previous_switch_positionA, $switch_positionA) = (undef, undef);
  my ($previous_switch_positionB, $switch_positionB) = (undef, undef);
  my @interXlengthsA = ();
  my @interXlengthsB = ();
  # print STDERR "progeny:  ", join(",  ", @$pgts2), "\n";
  # print STDERR "parent:   ", join(",  ", @$pgts1), "\n";
  while (my($i, $pgt1) = each @$pgts1) {
    my $par_pgt = $pgt1;	 # ->pgt();
    my $prog_pgt = $pgts2->[$i]; #->pgt();
    if ($par_pgt ne  'X'  and  $prog_pgt ne 'X') {
      # get the alleles of the two chromosomes of a homologous pair, here called A and B:
      my $progA = ($prog_pgt == 0 or $prog_pgt == 1)? 0 : 1; # 0 <-> prog A has ref allele.
      my $progB = ($prog_pgt == 0 or $prog_pgt == -1)? 0 : 1; # 0 <-> prog B has ref allele.
      if ($par_pgt == 0) {	# par: 0|0
	if ($prog_pgt == 0) {	# prog: 0|0
	  $Shom++;
	} elsif ($prog_pgt == 2) { # prog: 1|1
	  $Dhom++;
	}
      } elsif ($par_pgt == 1) {	# par 0|1

	$par_1count++;
	$par_het_count++;

	# relative to prog chrom A
	if ($progA == 0) {
	  $phaseA = 0;
	  #	  $A_0_count++;
	  if (defined $previous_phaseA  and  $phaseA != $previous_phaseA) { # recomb
	    $switch_countA++;
	    $switch_positionA = $par_het_count;
	    if (defined $previous_switch_positionA) {
	      my $intercross_length = $switch_positionA - $previous_switch_positionA;
	      push @interXlengthsA, $intercross_length;
	      $length_1_countA++ if(($switch_positionA - $previous_switch_positionA)  ==  1);
	    }
	    $previous_switch_positionA = $switch_positionA;
	  }
	  $previous_phaseA = $phaseA;
	} elsif ($progA == 1) {
	  $phaseA = 1;
	  #	  $A_1_count++;
	  if (defined $previous_phaseA  and  $phaseA != $previous_phaseA) {
	    $switch_countA++;
	    $switch_positionA = $par_het_count;
	    if (defined $previous_switch_positionA) {
	      my $intercross_length = $switch_positionA - $previous_switch_positionA;
	      push @interXlengthsA, $intercross_length;
	      # print STDERR  "0|1 A1 $switch_positionA ", $switch_positionA - $previous_switch_positionA, " $par_het_count\n";
	      $length_1_countA++ if($switch_positionA - $previous_switch_positionA  ==  1);
	    }
	    $previous_switch_positionA = $switch_positionA;
	  }
	  $previous_phaseA = $phaseA;
	}

	# relative to prog. chrom B.
	if ($progB == 0) {
	  $phaseB = 0;
	  #	  $B_0_count++;
	  if (defined $previous_phaseB  and  $phaseB != $previous_phaseB) {
	    $switch_countB++;
	    $switch_positionB = $par_het_count;
	    if (defined $previous_switch_positionB) {
	      my $intercross_length = $switch_positionB - $previous_switch_positionB;
	      push @interXlengthsB, $intercross_length;
	      $length_1_countB++ if($switch_positionB - $previous_switch_positionB  ==  1);
	    }
	    $previous_switch_positionB = $switch_positionB;
	    
	  }
	  $previous_phaseB = $phaseB;
	} elsif ($progB == 1) {
	  $phaseB = 1;

	  #	  $B_1_count++;
	  if (defined $previous_phaseB  and  $phaseB != $previous_phaseB) {
	    $switch_countB++;
	    $switch_positionB = $par_het_count;
	    if (defined $previous_switch_positionB) {
	      my $intercross_length = $switch_positionB - $previous_switch_positionB;
	      push @interXlengthsB, $intercross_length;
	      $length_1_countB++ if($switch_positionB - $previous_switch_positionB  ==  1);
	    }
	    $previous_switch_positionB = $switch_positionB;
	  }
	  $previous_phaseB = $phaseB;
	}

      } elsif ($par_pgt == -1) {	# par 1|0
	$par_2count++;
	$par_het_count++;

	# relative to prog. chrom. A.
	if ($progA == 0) {
	  $phaseA = 1;
	  if (defined $previous_phaseA  and  $phaseA != $previous_phaseA) {
	    $switch_countA++;
	    $switch_positionA = $par_het_count;
	    if (defined $previous_switch_positionA) {
	      my $intercross_length = $switch_positionA - $previous_switch_positionA;
	      push @interXlengthsA, $intercross_length;
	      $length_1_countA++ if($switch_positionA - $previous_switch_positionA  ==  1);
	    }
	    $previous_switch_positionA = $switch_positionA;
	  }
	  $previous_phaseA = $phaseA;
	} elsif ($progA == 1) {
	  $phaseA = 0;
	  if (defined $previous_phaseA  and  $phaseA != $previous_phaseA) {
	    $switch_countA++;
	    $switch_positionA = $par_het_count;
	    if (defined $previous_switch_positionA) {
	      my $intercross_length = $switch_positionA - $previous_switch_positionA;
	      push @interXlengthsA, $intercross_length;
	      $length_1_countA++ if($switch_positionA - $previous_switch_positionA  ==  1);
	    }
	    $previous_switch_positionA = $switch_positionA;
	  }
	  $previous_phaseA = $phaseA;
	}
	if ($progB == 0) {
	  $phaseB = 1;
	  if (defined $previous_phaseB  and  $phaseB != $previous_phaseB) {
	    $switch_countB++;
	    $switch_positionB = $par_het_count;
	    if (defined $previous_switch_positionB) {
	      my $intercross_length = $switch_positionB - $previous_switch_positionB;
	      push @interXlengthsB, $intercross_length;
	      $length_1_countB++ if($switch_positionB - $previous_switch_positionB  ==  1);
	    }
	    $previous_switch_positionB = $switch_positionB;
	  }
	  $previous_phaseB = $phaseB;
	} elsif ($progB == 1) {
	  $phaseB = 0;
	  if (defined $previous_phaseB  and  $phaseB != $previous_phaseB) {
	    $switch_countB++;
	    $switch_positionB = $par_het_count;
	    if (defined $previous_switch_positionB) {
	      my $intercross_length = $switch_positionB - $previous_switch_positionB;
	      push @interXlengthsB, $intercross_length;
	      $length_1_countB++ if($switch_positionB - $previous_switch_positionB  ==  1);
	    }
	    $previous_switch_positionB = $switch_positionB;
	  }
	  $previous_phaseB = $phaseB;
	}
      } elsif ( $par_pgt == 2) { # par: 1|1
	if ($prog_pgt == 0) {	 # prog: 0|0
	  $Dhom++;
	} elsif ($prog_pgt == 2) { # prog: 1|1
	  $Shom++;
	}
      }
    }
    # print STDERR "$prog_pgt  $par_pgt  $Dhom $Shom\n"
  }
  #  print "A  ", join(" ", @interXlengthsA), "\n";
  my $n_dub_A =  count_dubious_Xs(\@interXlengthsA) . "\n";
  # print "B  ", join(" ", @interXlengthsB), "\n";
  my $n_dub_B = count_dubious_Xs(\@interXlengthsB) . "\n";
    
  return ($Dhom, $Shom, $switch_countA, $switch_countB, $par_het_count, $n_dub_A, \@interXlengthsA, $n_dub_B, \@interXlengthsB, );
}

sub print_help{
  my $stream = shift;
  print STDOUT "usage:  phased_parents  -pdsgs <pdsgs_file>  -ped <pedigree_file>  [options]\n";
  print STDOUT "-pdsgs_file   <pdsgs file>  (pdsgs file must have phased dosages, i.e. 1, -1 for heterozygs, and chromosomes on line 2)\n";
  print STDOUT "-pedigree_file  <pedigree file> (Accession, Female parent, Male parent  ids in columns 1, 2, and 3)\n";
  print STDOUT "-output_file <output filename> (default: 'pp.out')\n";
  print STDOUT "-n_rand_parents  N  (do N randomly chosen 'parents' per accession in addition to pedigree parents (default: 0)\n";
  print STDOUT "-reverse  (causes calculation of reversed parent-offspring relationship,\n             i.e. P as offspring of A, as well as usual A as offspring of P (default: false)\n";
  print STDOUT "-help  Print this message.\n";
}


  sub count_L1s{
    my $interXls = shift;
    my $count = 0;
    for my $iXl (@$interXls) {
      $count++ if($iXl == 1);
    }
    return $count;
  }

sub count_dubious_Xs{ # the Xs which could be eliminated if flip single isolated phases (not allowed to flip two or more in a row)
  my $interXls = shift;
  my $count = 0;
  my $n = scalar @$interXls;
  my $prev_iXl = undef;
  my $ones_in_a_row = 0;
  my $n_dub;
  for (my $i=0; $i<$n; $i++) {
    my $iXl = $interXls->[$i];
    if ($iXl == 1) {
      #if (defined $prev_iXl  and  $prev_iXl == 1) {
      $ones_in_a_row++;
      #}
    } else {
      if (defined $prev_iXl  and  $prev_iXl == 1) { # previous one was the last in a series of (1 or more) 1's

	$n_dub = $ones_in_a_row + ($ones_in_a_row % 2);
	$count +=  $n_dub;

	#print "$prev_iXl  $iXl  $ones_in_a_row   $n_dub  $count\n";
      }
      $ones_in_a_row = 0;
    }
    $prev_iXl = $iXl;
  }
  if (defined $prev_iXl  and  $prev_iXl == 1) {
    $n_dub = $ones_in_a_row + ($ones_in_a_row % 2);
    $count +=  $n_dub;
    #print "XX $prev_iXl  ---  $ones_in_a_row  $n_dub  $count\n";
  }
  return $count;
}
