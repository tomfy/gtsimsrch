#!/usr/bin/perl -w
use strict;

my %id_p = ();
my ($nNtotal, $nQtotal) =  (0, 0); # (499000, 1000); # (99800, 200); #
# $nNtotal *= 4; # ???
# $nQtotal *= 4; # ???
my $count = 0;
my $bonf_threshold = 1e-5;
my ($nsigN, $nsigQ) = (0, 0);
while(<>){ # read an .mlma gwas output file with sp id in col 2, p in col 9
  next if(/^\s*(#|Chr)/);
  my @cols = split(" ", $_);
  my ($id, $p) = @cols[1,8];
  $nQtotal++ if($id =~ /Q/);
  $nNtotal++ if($id =~ /Null/);
  $id .= '_' . $count;
  $id_p{$id} = $p;
  $count++;
}

my @sorted_ids = sort { $id_p{$b} <=> $id_p{$a} } keys %id_p; # sort by value (p) decreasing.

my ($nN, $nQ) = (0, 0); # count the number of Null and Q snps with $p greater than various thresholds
# i.e. $nN is false positives, $nQ is true positives
# 
#for my $id (@sorted_ids){
  while(my($i, $id) = each @sorted_ids){
   # my $next_id = ($i < $#sorted_ids)? $sorted_ids[$i+1] : $id;
    my $ppp =  $id_p{$id}; # 0.5*($id_p{$id} + $id_p{$next_id});
  #  print STDERR "id: $id\n";
   if($id =~ /Null/){
     $nN++;
     $nsigN++ if($id_p{$id} < $bonf_threshold);
  }elsif($id =~ /Q/){
    $nQ++;
    $nsigQ++ if($id_p{$id} < $bonf_threshold);
  }else{
    print STDERR "unexpected id encountered: $id \n";
    die;
  }
   my $fnr = ($nQ)/$nQtotal;
   my $fpr = ($nNtotal -$nN)/$nNtotal;
  print "$nN  $nQ ", $nNtotal-$nN, "  $fpr  $fnr   $id  $ppp\n";
 }
print "# nsigN: $nsigN  nsigQ  $nsigQ \n";
