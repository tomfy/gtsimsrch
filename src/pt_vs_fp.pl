#!/usr/bin/perl -w
use strict;
use List::Util qw(min max sum);

# compare the analyses by pedigree_test and find_parents / best_triples

my $pedigree_test_file = shift;

my $best_triples_file = shift;

open my $fhpt, "<", "$pedigree_test_file";

# my @pt_cats = ({}, {}, {}, {}, {});
my %ptids_cat = ();
my $cFr = 0.05;
my $cHgmr = 0.04;
my $cD = 0.03;

while(<$fhpt>){
  next if(/^\s*#/);
  my @cols = split(" ", $_);
  my ($accid, $Fid, $Mid, $Fhgmr, $Fr, $Mhgmr, $Mr, $d) = @cols[0, 2,3, 7,9,11,13, 17];
  my $the_cat; # 1,2,3,4,5
  if($Fid eq 'NA'  or  $Mid eq 'NA'){

  }else{
  if($Fid eq $Mid){ # self according to pedigree file
    if($Fr <= $cFr  and  $d <= $cD){
      $the_cat = 1;
    }else{
      $the_cat = 2;
    }
  }else{ # bip according to pedigree file
    if($Fr > $cFr  and  $d <= $cD){
      $the_cat = 3;
    }elsif($Fr > $cFr  and  $Fhgmr < $cHgmr  and $Mhgmr < $cHgmr){
      $the_cat = 4;
    }else{
      $the_cat = 5;
    }
  }
}
  my $parents_pair = ($Fid le $Mid)? "$Fid $Mid" : "$Mid $Fid";
   $ptids_cat{"$accid $parents_pair"} = $the_cat;
  #$pt_cats[$the_cat-1]->{"$accid $parents_pair"} = 1;
}
close $fhpt;



open my $fhbt, "<", "$best_triples_file";

#my $bt_cats = [{}, {}, {}, {}];
my %btids_cat = ();
my $cDa = 0.02;
my $cDb = 0.09;

while(<$fhbt>){
  next if(/^\s*#/);
  my @cols = split(" ", $_);
  my ($accid, $Fid, $Mid, $Fhgmr, $Fr, $Mhgmr, $Mr, $d1) = @cols[0, 3,4, 9,11,13,15, 19];
 my $the_cat; # 1,2,3,4,5
  if(scalar @cols > 43){
    my $d2 = $cols[42];
    
 
  if($d1 > $cD){
    $the_cat = 1;
  }elsif($d2 > $cDb){
    $the_cat = 2;
  }elsif($d2 > $cDa  and  (2*$d1 <= $d2)){
      $the_cat = 3;
    }else{
      $the_cat = 4;
    }
  }else{ # fp found only one solution
    $the_cat = 5;
  }
 # }
  my $parents_pair = ($Fid le $Mid)? "$Fid $Mid" : "$Mid $Fid";
  $btids_cat{"$accid $parents_pair"} = $the_cat;
  # $bt_cats[$the_cat-1]->{"$accid $parents_pair"} = 1;
  
}
close $fhbt;

my @table = ([0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0], [0,0,0,0,0,0]);

my $bt_not_pt_count = 0;
while(my ($btids, $btcat) = each %btids_cat){
#  print STDERR "$btids  $btcat\n";
  if(exists $ptids_cat{$btids}){
    my $ptcat = $ptids_cat{$btids};
  #  print STDERR "    $ptcat \n";
    $table[$ptcat]->[$btcat]++;
  }else{
    $table[0]->[$btcat]++; # absent from pedigree_test output
  }
}

while(my ($ptids, $ptcat) = each %ptids_cat){

 # print STDERR "xxx: $ptids  $ptcat \n";
  if(!exists $btids_cat{$ptids}){
    $table[$ptcat]->[0]++;
  }
}
my ($table_sum_a, $table_sum_b) = (0, 0);
my @col_sums = (0, 0, 0, 0, 0, 0);
while(my($pcat, $counts) = each @table){
  printf STDOUT "ptcat: %2d    ", $pcat;
  while(my($bcat, $acount) = each @$counts){
    printf STDOUT "%6d  ", $acount;
  $col_sums[$bcat] += $acount;
}
  my $ptcat_sum = sum(@$counts);
  $table_sum_a += $ptcat_sum;
  printf STDOUT "   row sum:  %6d", $ptcat_sum;
  print "\n";
}

printf STDOUT "colsums:     ";
for my $acolsum (@col_sums){
  printf STDOUT "%6d  ", $acolsum;
}
$table_sum_b = sum(@col_sums);
if($table_sum_a == $table_sum_b){
  printf STDOUT " table sum:  %6d", $table_sum_a;
}else{
  printf "table sums not equal?  $table_sum_a  $table_sum_b \n";
}
print "\n";

  


