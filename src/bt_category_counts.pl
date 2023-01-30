#!/usr/bin/perl -w
use strict;

my $d1_max = shift // 0.03;
my $d2_1 = shift // 0.09;
my $d2_2 = shift // 0.018;
my $dratio = shift // 2;
my $d1_col = shift // 20;
my $d2_col = shift // 43;
$d1_col--; $d2_col--; # now zero-based

my @ok1s = ();
my @ok2s = ();
my @ok3s = ();
my @ok4s = ();
my @ambigs = ();
my @bads = ();

my ($n1, $n2, $n3, $n4, $n5, $n6) = (0, 0, 0, 0, 0, 0);
while(my $line = <>){
  my @cols = split(" ", $line);
  my $str = $cols[0] . " ";
  $str .= join(" ", @cols[2,3,4]) . " " . $cols[19];

  my $d1 = $cols[$d1_col];
  if($cols[5] > 1){ # more than 1 soln
    my $d2 = $cols[$d2_col];
      $str .= " " . join(" ", @cols[25,26,27]) . " " . $cols[42];
  if($d1 >= $d1_max){ # d1 is bad
    $n1++;
    push @bads, $str;
  }elsif($dratio*$d1 > $d2){ # d2 almost as good as d1 (within factor of 2)
    $n5++;
    push @ambigs, $str;
  }elsif($d2 > $d2_1){ # d1 good, d2 very bad
    $n2++;
    push @ok2s, $str;
  }elsif($d2 > $d2_2){ # d1 good, d2 bad
    $n3++;
    push @ok3s, $str;
  }else{ # d1 good, d2 also good, but at least twice as large, so significantly less good.
    $n4++;
    push @ok4s, $str;
  }
}else{ # only one solution found
  if($d1 >= $d1_max){ # d1 bad (so zero good solutions found)
    $n1++;
    push @bads, $str;
  }else{
    $n6++; # d1 good, no other solutions
    push @ok1s, $str;
  }
}
}
my $total = $n1 + $n2 + $n3 + $n4 + $n5 + $n6;

#print "$total  $n1  $n2  $n3  $n4  $n5  $n6\n";
print "d1max: $d1_max  d2_1: $d2_1  d2_2: $d2_2  dratio: $dratio \n";
print "  total   ok1    ok2    ok3    ok4  ambig    bad \n";
printf(STDOUT "%6i %6i %6i %6i %6i %6i %6i \n",
       $total, $n6, $n2, $n3, $n4,  $n5, $n1);

open my $fh, ">", "ok1.out";
for (@ok1s){
  print $fh "$_\n";
}
close $fh;

open $fh, ">", "ok2.out";
for (@ok2s){
  print $fh "$_\n";
}
close $fh;

open $fh, ">", "ok3.out";
for (@ok3s){
  print $fh "$_\n";
}
close $fh;

open $fh, ">", "ok4.out";
for (@ok4s){
  print $fh "$_\n";
}
close $fh;

open $fh, ">", "ambig.out";
for (@ambigs){
  print $fh "$_\n";
}
close $fh;

open $fh, ">", "bad.out";
for (@bads){
  print $fh "$_\n";
}
close $fh;
