#!/usr/bin/perl -w
use strict;
use List::Util qw(min max sum shuffle);

# read in a vcf file
# if it has phased genotypes, analyze to find
# parent-progeny pairs.

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
print STDERR "# libdir: $libdir \n";

use Chromosome;
#use Marker;
my $do_reverse = 0;

my $pedigree_file = shift // undef;
my $rand_parents_to_do = shift // 4;

##########################################################
# read pedigree file and store accession-parent pairs ####
##########################################################
my %A_Fpar = ();
my %A_Mpar = ();
if (defined $pedigree_file) {
  open my $fhped, "<", "$pedigree_file";
  while (<$fhped>) {
    my @cols = split(" ", $_);
    my ($A, $Fpar, $Mpar) = @cols[0,1,2]; # [-3,-2,-1];
    #print STDERR "$A  $Fpar $Mpar\n";
    $A_Fpar{$A} = $Fpar if($Fpar ne 'NA');
    $A_Mpar{$A} = $Mpar if($Mpar ne 'NA');

  }
  close $fhped;
}
print STDERR  "# Pedigree information from file  $pedigree_file  stored.\n";
print STDERR "# ", scalar keys %A_Fpar, " female parents;  ", scalar keys %A_Mpar, " male parents.\n";
#sleep(2);
##########################################################
# Done reading pedigree file ######################## ####
##########################################################


##############################################
# read and store vcf file ####################
##############################################
my $max_markers = 5000000;
my @acc_ids = ();
my %acc_chroms = (); # keys: acc ids, values: arrayref of chromosome object;
my %chrom_numbers = ();

while (my $line = <>) {		# read accession ids from vcf file
  next if($line =~ /^\s*##/);
  if ($line =~ /^\s*#/) {
    @acc_ids = split(" ", $line);
    @acc_ids = @acc_ids[9..$#acc_ids]; #
    print STDERR "# number of accession ids in vcf file: ", scalar @acc_ids, "\n";
    last;
  }
}
# print 'YYY ',  scalar @acc_ids, "\n";
for my $acc_id (@acc_ids) {
  $acc_chroms{$acc_id} = [];
}

# read the vcf file and store the phased genotype information
my $marker_count = 0;
while (my $line = <>) # each pass through this loop processes one marker, all accessions
  {
    my @marker_gts = split(" ", $line); # these are the (phased) genotypes for the various accessions, this marker
    # things like  0|0:0:1,0,0
    # i.e. colon separated fields. here we require the first of these to be '0|0' '0|1' '1|0', '1|1'
    my $i_chrom = shift @marker_gts;
    $chrom_numbers{$i_chrom} = 1;
    my $position = shift @marker_gts;
    my $mrkr_id = shift @marker_gts;
    my $ref = shift @marker_gts;
    my $alt = shift @marker_gts;
    my $qual = shift @marker_gts;
    my $filter = shift @marker_gts;
    my $info = shift @marker_gts;
    my $format = shift @marker_gts;
    # print STDERR  scalar @acc_ids, "  ", scalar @marker_gts, "\n";
    while (my($i, $accid) = each @acc_ids) { # loop over the accessions
      # print STDERR "i: $i  accid:  $accid\n";
      if (!exists $acc_chroms{$accid}->[$i_chrom]) { # check whether Chromosome obj. exists for this accession and chromosome number
	$acc_chroms{$accid}->[$i_chrom] = Chromosome->new({i_chrom => $i_chrom, genotypes => []});
      }
      my $the_chrom = $acc_chroms{$accid}->[$i_chrom];
      # print ref $the_chrom, "\n";
      # print $the_chrom->i_chrom(), "\n";
      # print scalar @{$the_chrom->genotypes()}, "\n";

      my $the_pgt = '4';	# 4 means missing data
      my $gt_info = $marker_gts[$i];
      my @fields = split(':', $gt_info);
      my @gps = split(',', $fields[2]); # est probabilities for different 
      if (max(@gps) < 0.9) {
	# print STDERR join(':', @fields), "\n";
	$the_pgt = 4;		# low quality gt, set as missing.
      } else {
	if ($fields[0] eq '0|0') {
	  $the_pgt = 0;
	} elsif ($fields[0] eq '0|1') {
	  $the_pgt = 1;
	} elsif ($fields[0] eq '1|0') {
	  $the_pgt = 2;
	} elsif ($fields[0] eq '1|1') {
	  $the_pgt = 3;
	}
      }
      #    my $the_marker = Marker->new({ pgt => $the_pgt});
      #    $the_chrom->add_marker($the_marker);
      $the_chrom->add_genotype($the_pgt);
      # print STDERR "Added genotypes for $accid, chrom $i_chrom.\n";
    }
    $marker_count++;
    printf STDERR "Markers read from vcf file so far: $marker_count \n" if($marker_count % 100 == 0);
    last if($marker_count >= $max_markers);
  }
printf STDERR "Done storing markers. Markers stored:  $marker_count \n";
###########################################
# done storing info in vcf file. ##########
###########################################

# for my $anid (@acc_ids){
#   print STDERR "$anid   ", $A_Fpar{$anid} // 'X', "  ", $A_Mpar{$anid} // 'Y', "\n";
# }

my @chroms = sort {$a <=> $b} keys %chrom_numbers;
my $n_chrom_pairs_analyzed_forward = 0;
# for my $i_chrom (@chroms) { # loop over chromosomes
while (my ($i, $progid) = each @acc_ids) { # loop over accessions considered as progeny
  my $ped_Fpar = $A_Fpar{$progid} // 'unknown';
  my $ped_Mpar = $A_Mpar{$progid} // 'unknown';
  my %cand_parent_ids = ();
  if ($ped_Fpar ne 'unknown') {
    $cand_parent_ids{$ped_Fpar} = 'F';
  }
  if ($ped_Mpar ne 'unknown') {
    $cand_parent_ids{$ped_Mpar} = 'M';
  }
  my $rand_parents_added = 0;
  my @shuffled_ids = shuffle(@acc_ids);
  for my $anid (@shuffled_ids) { # add some other randomly chosen parents
    if (exists $acc_chroms{$anid}) {
      if ($anid ne $progid  and  $anid ne $ped_Fpar  and  $anid ne $ped_Mpar) {
	$cand_parent_ids{$anid} = 'R';
	$rand_parents_added++;
      }
    }
    last if($rand_parents_added >= $rand_parents_to_do);
  }
  for my $i_chrom (@chroms) {	# loop over chromosomes
    my $pgts1 = $acc_chroms{$progid}->[$i_chrom]->genotypes();
    if ($ped_Fpar eq 'unknown'  and  $ped_Mpar eq 'unknown') {
      print STDERR "Accession $progid has both parent unknown in pedigree file.\n";
      next;
    }

    while (my($parid, $type) = each %cand_parent_ids) {

      next if($parid eq $progid);
      #my $pFM = $type;		# $pF . $pM;
      next if(!exists $acc_chroms{$parid});
      # if ($pFM ne 'XX') {
      #print STDERR "$i_chrom  $progid  $pFM   $ped_Fpar  $ped_Mpar \n";
      # forward direction - $parid is parent
      my $pgts2 = $acc_chroms{$parid}->[$i_chrom]->genotypes();
      my ($Do, $So, $XA, $XB, $parent_het_count, $length1count_A, $length1count_B) = analyze_pgts_pair($pgts2, $pgts1); # $pgts1: parent, $pgts2: progeny
      $n_chrom_pairs_analyzed_forward++;
      print STDERR "Chromosome parent-progeny pairs analyzed: $n_chrom_pairs_analyzed_forward \n" if($n_chrom_pairs_analyzed_forward % 100 == 0);
      my $hgmr_denom = $Do + $So;
      # $par_01_count is number of 0|1 genotypes in parent ($progid)
      # $par_10_count is number of 1|0 genotypes in parent ($progid)
      print "$progid  $parid  $i_chrom  ", ($hgmr_denom > 0)? $Do/$hgmr_denom : '-1', "   $XA $XB   ";
      if ($XA < $XB) {
	print "$XA $XB  ";	# $length1count_A $length1count_B  ";
      } else {
	print "$XB $XA  ";	#  $length1count_B $length1count_A  ";
      }		       # , min($XA, $XB), "  ", max($XA, $XB), "   ", 
      #      "   $par_01_count $par_10_count ",
      print " $parent_het_count  $type  forward\n";

      if ($do_reverse) { # now with $progid as parent, $parid as progeny
	my ($Do_rev, $So_rev, $XA_rev, $XB_rev, $parent_het_count_rev, $L1count_A_rev, $L1count_B_rev) = analyze_pgts_pair($pgts1, $pgts2);
	my $hgmr_rev = ($So_rev > 0)? $Do_rev/($Do_rev + $So_rev) : -1;
	print "$parid  $progid  $i_chrom  $hgmr_rev  $XA_rev $XB_rev   ";
	if ($XA_rev < $XB_rev) {
	  print "$XA_rev $XB_rev  "; # $length1count_A $length1count_B  ";
	} else {
	  print "$XB_rev $XA_rev  "; #  $length1count_B $length1count_A  ";
	}	       # , min($XA, $XB), "  ", max($XA, $XB), "   ", 
	#      "   $par_01_count $par_10_count ",
	print " $parent_het_count_rev reverse\n";
      }
      # }
    }
  }
}

#  ****************************************************************************

sub analyze_pgts_pair{ # consider $pgts1 as parent, $pgts2 as progeny.
  my $pgts1 = shift;
  my $pgts2 = shift;
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
  #  my ($A_0_count, $A_1_count, $B_0_count, $B_1_count) = (0, 0, 0, 0);
  # print join('', map($_->pgt(), @$pgts1)), "\n";
  # print join('', map($_->pgt(), @$pgts2)), "\n";
  my ($previous_switch_positionA, $switch_positionA) = (undef, undef);
  my ($previous_switch_positionB, $switch_positionB) = (undef, undef);
  while (my($i, $pgt1) = each @$pgts1) {
    my $par_pgt = $pgt1;	 # ->pgt();
    my $prog_pgt = $pgts2->[$i]; #->pgt();
    # print "$par_pgt  $prog_pgt \n";
    if ($par_pgt != 4  and  $prog_pgt != 4) {
      my $progA = ($prog_pgt == 0 or $prog_pgt == 1)? 0 : 1; # 0 <-> prog A has 0 (ref)
      my $progB = ($prog_pgt == 0 or $prog_pgt == 2)? 0 : 1; # 0 <-> prog B has 0 (ref)
      if ($par_pgt == 0) {	# par: 0|0
	if ($prog_pgt == 0) {	# prog: 0|0
	  $Shom++;
	} elsif ($prog_pgt == 3) { # prog: 1|1
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
	      # print STDERR  "0|1 A0 $switch_positionA  ", $switch_positionA - $previous_switch_positionA, " $par_het_count\n";
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
	      # print STDERR  "0|1 B0 $switch_positionB ", $switch_positionB - $previous_switch_positionB, "  $par_het_count\n";
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
	      # print STDERR  "0|1 B1 $switch_positionB ", $switch_positionB - $previous_switch_positionB, " $par_het_count\n";
	      $length_1_countB++ if($switch_positionB - $previous_switch_positionB  ==  1);
	    }
	    $previous_switch_positionB = $switch_positionB;
	  }
	  $previous_phaseB = $phaseB;
	}

	#			print "X $par_pgt  $prog_pgt   $progA $progB  $previous_phaseA $previous_phaseB   $switch_countA $switch_countB \n";

      } elsif ($par_pgt == 2) {	# par 1|0
	$par_2count++;
	$par_het_count++;

	# relative to prog. chrom. A.
	if ($progA == 0) {
	  $phaseA = 1;
	  if (defined $previous_phaseA  and  $phaseA != $previous_phaseA) {
	    $switch_countA++;
	    $switch_positionA = $par_het_count;
	    if (defined $previous_switch_positionA) {
	      # print STDERR  "1|0 A0 $switch_positionA ", $switch_positionA - $previous_switch_positionA, " $par_het_count\n";
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
	      # print STDERR  "1|0 A1 $switch_positionA ", $switch_positionA - $previous_switch_positionA, " $par_het_count\n";
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
	      # print STDERR  "1|0 B0 $switch_positionB ", $switch_positionB - $previous_switch_positionB, " $par_het_count\n";
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
	      # print STDERR  "1|0 B1 $switch_positionB ", $switch_positionB - $previous_switch_positionB, " $par_het_count\n";
	      $length_1_countB++ if($switch_positionB - $previous_switch_positionB  ==  1);
	    }
	    $previous_switch_positionB = $switch_positionB;
	  }
	  $previous_phaseB = $phaseB;
	}
      } elsif ( $par_pgt == 3) { # par: 1|1
	if ($prog_pgt == 0) {	 # prog: 0|0
	  $Dhom++;
	} elsif ($prog_pgt == 3) { # prog: 1|1
	  $Shom++;
	}
      }
      #  print STDERR "$par_pgt  $progA $progB  " , $phaseA // '-', "  $switch_countA  ", $phaseB // '-', "  $switch_countB\n";
    }
  }
  return ($Dhom, $Shom, $switch_countA, $switch_countB, $par_het_count, $length_1_countA, $length_1_countB);
}
