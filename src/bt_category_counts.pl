#!/usr/bin/perl -w
use strict;

my $d1_max = shift // 0.03;
my $d2_1 = shift // 0.09;
my $d2_2 = shift // 0.018;
my $dratio = shift // 2;
my $d1_col = shift // 20;
my $d2_col = shift // 43;
$d1_col--; $d2_col--; # now zero-based

my ($n1, $n2, $n3, $n4, $n5, $n6) = (0, 0, 0, 0, 0, 0);
while(my $line = <>){
  my @cols = split(" ", $line); 
  my $d1 = $cols[$d1_col];
  if($cols[5] > 1){ # more than 1 soln
  my $d2 = $cols[$d2_col];
  if($d1 >= $d1_max){ # d1 is bad
    $n1++;
  }elsif($dratio*$d1 > $d2){ # d2 almost as good as d1 (within factor of 2)
    $n5++;
  }elsif($d2 > $d2_1){ # d1 good, d2 very bad
    $n2++;
  }elsif($d2 > $d2_2){ # d1 good, d2 bad
    $n3++;
  }else{ # d1 good, d2 also good, but at least twice as large.
    $n4++;
  }
}else{ # only one solution found
  if($d1 >= $d1_max){ # d1 bad (so zero good solutions found)
    $n1++;
  }else{
    $n6++; # d1 good, no other solutions
  }
}
}
my $total = $n1 + $n2 + $n3 + $n4 + $n5 + $n6;

#print "$total  $n1  $n2  $n3  $n4  $n5  $n6\n";
print "d1max: $d1_max  d2_1: $d2_1  d2_2: $d2_2  dratio: $dratio \n";
print "  total   ok1    ok2    ok3    ok4  ambig    bad \n";
printf(STDOUT "%6i %6i %6i %6i %6i %6i %6i \n", $total, $n6, $n2, $n3, $n4,  $n5, $n1);
