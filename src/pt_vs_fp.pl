#!/usr/bin/perl -w
use strict;
use List::Util qw(min max sum);

# compare the analyses by pedigree_test and find_parents / best_triples

my $pedigree_test_file = shift;

my $best_triples_file = shift;

my $n_pt_categories = 6;

open my $fhpt, "<", "$pedigree_test_file";

# my @pt_cats = ({}, {}, {}, {}, {});
my %ped_lines = ();
my %ptids_cat = ();
my $cFr = 0.05;
my $cHgmr = 0.04;
my $cD = 0.03;

my %ptaccid = ();
my %btaccid = ();

while (<$fhpt>) {
  next if(/^\s*#/);
  my @cols = split(" ", $_);
  my ($accid, $Fid, $Mid, $Fhgmr, $Fr, $Mhgmr, $Mr, $d) = @cols[0, 2,3, 7,9,11,13, 17];
  
  $ptaccid{$accid} = 1;
  my $the_cat;				 # 1,2,3,4,5
  if ($Fid ne 'NA'  and  $Mid eq 'NA') { # F parent only: good hgmr (large/small R), bad hgmr
    # if ($Fhgmr <= $cHgmr) {
    #   if ($Fr <= $cFr) {
    # 	$the_cat = 7;		# F only, probable self
    #   } else {
    # 	$the_cat = 6;		# F only, F ok, M unknown
    #   }
    # } else {
    #   $the_cat = 8;		# F only, probably wrong
    # }
    $the_cat = 0;
  } elsif ($Fid eq 'NA' and $Mid ne 'NA') { # M parent only: 
    # if ($Mhgmr <= $cHgmr) {
    #   if ($Mhgmr <= $cFr) {
    # 	$the_cat = 10;		# M only, probable self
    #   } else {
    # 	$the_cat = 9;		# M only, F ok, M unknown
    #   }
    # } else {
    #   $the_cat = 11;		# M only, probably wrong
    # }
    $the_cat = 0;
  } elsif ($Fid eq 'NA' and $Mid eq 'NA') {
    $the_cat = 0;
  } else {
    if ($Fid eq $Mid) {		# self according to pedigree file
      if ($Fr <= $cFr  and  $d <= $cD) {
	$the_cat = 1;
      } else {
	$the_cat = 2;
      }
    } else {			# bip according to pedigree file
      if ($Fr > $cFr  and  $d <= $cD) {
	$the_cat = 3;
      } elsif ($Fr > $cFr  and  $Fhgmr < $cHgmr  and $Mhgmr < $cHgmr) {
	$the_cat = 4;
      } else {
	$the_cat = 5;
      }
    }
  }
  my $parents_pair = ($Fid le $Mid)? "$Fid $Mid" : "$Mid $Fid";
  $ptids_cat{"$accid $parents_pair"} = $the_cat;
  $ped_lines{"$accid $parents_pair"} = $_;
  #$pt_cats[$the_cat-1]->{"$accid $parents_pair"} = 1;
}
close $fhpt;



open my $fhbt, "<", "$best_triples_file";

#my $bt_cats = [{}, {}, {}, {}];
my %btids_cat = ();
my $cDa = 0.02;
my $cDb = 0.09;

while (<$fhbt>) {
  next if(/^\s*#/);
  my @cols = split(" ", $_);
  my ($accid, $Fid, $Mid, $Fhgmr, $Fr, $Mhgmr, $Mr, $d1) = @cols[0, 3,4, 9,11,13,15, 19];
 
  $btaccid{$accid} = 1;
  my $the_cat;			# 1,2,3,4,5
  if (scalar @cols > 43) {
    my $d2 = $cols[42];
   
 
    if ($d1 > $cD) {
      $the_cat = 1;
    } elsif ($d2 > $cDb) {
      $the_cat = 2;
    } elsif ($d2 > $cDa  and  (2*$d1 <= $d2)) {
      $the_cat = 3;
    } else {
      $the_cat = 4;
    }

    my ($Fid2, $Mid2) = @cols[26,27];
    my $parents_pair2 = ($Fid2 le $Mid2)? "$Fid2 $Mid2" : "$Mid2 $Fid2";
    if ($d2 <= $cD) {
        $btids_cat{"$accid $parents_pair2"} = 7;
    } else {
        $btids_cat{"$accid $parents_pair2"} = 8;
    }
  } else { # fp found only one solution, i.e. only one cand. parent, so only solution is a self.
    if ($d1 <= $cD) {
      $the_cat = 5;
    } else {
      $the_cat = 6;
    }
  }
  # }
  my $parents_pair = ($Fid le $Mid)? "$Fid $Mid" : "$Mid $Fid";
  $btids_cat{"$accid $parents_pair"} = $the_cat;

  
  # $bt_cats[$the_cat-1]->{"$accid $parents_pair"} = 1;
  
}
close $fhbt;

my @table = ();
#	     [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0],
#	     [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0]
#	    );
for (0..$n_pt_categories-1) {
  my @zeros = ((0) x 9);
  push @table, \@zeros;
}
my $bt_not_pt_count = 0;
while (my ($btids, $btcat) = each %btids_cat) {
  #  print STDERR "$btids  $btcat\n";
  if (exists $ptids_cat{$btids}) {
    my $ptcat = $ptids_cat{$btids};
    #  print STDERR "    $ptcat \n";
    if($btcat == 6  and $ptcat == 2){
      #print STDERR "$btids\n";
      print STDERR $ped_lines{$btids};
    }
    $table[$ptcat]->[$btcat]++;
  } else {
    $table[0]->[$btcat]++;	# absent from pedigree_test output
  }
}

while (my ($ptids, $ptcat) = each %ptids_cat) {

  # print STDERR "xxx: $ptids  $ptcat \n";
  if (!exists $btids_cat{$ptids}) {
    $table[$ptcat]->[0]++;
  }
}
my ($table_sum_a, $table_sum_b) = (0, 0);
my @col_sums = ((0) x 9);
while (my($pcat, $counts) = each @table) {
  printf STDOUT "ptcat: %2d    ", $pcat;
  while (my($bcat, $acount) = each @$counts) {
    printf STDOUT "%6d  ", $acount;
    $col_sums[$bcat] += $acount;
  }
  my $ptcat_sum = sum(@$counts);
  $table_sum_a += $ptcat_sum;
  printf STDOUT "   row sum:  %6d", $ptcat_sum;
  print "\n";
}

printf STDOUT "colsums:     ";
for my $acolsum (@col_sums) {
  printf STDOUT "%6d  ", $acolsum;
}
$table_sum_b = sum(@col_sums);
if ($table_sum_a == $table_sum_b) {
  printf STDOUT " table sum:  %6d", $table_sum_a;
} else {
  printf "table sums not equal?  $table_sum_a  $table_sum_b \n";
}
print "\n";

my $pt_only_count = 0;
my $both_count = 0;
for my $ptaid (keys %ptaccid) {
  if (exists $btaccid{$ptaid}) {
    $both_count++;
  } else {
    $pt_only_count++;
  }
}
my $pt_count = scalar keys %ptaccid;
my $fp_count = scalar keys %btaccid;
my $fp_only_count = $fp_count - $both_count;
print "pedigree_test count:              $pt_count\n";
print "find_parents/best_triple count:   $fp_count\n";
print "counts, ptonly, both, fponly: $pt_only_count  $both_count  $fp_only_count \n"; 
  


