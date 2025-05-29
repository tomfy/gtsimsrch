#!/usr/bin/perl -w
use strict;
use List::Util qw(min max sum);

my $OK = 'OK';
my $bad_colon_count = 'bad colon count';
my $noGT = 'bad GT length';
my $noDS = 'bad DS length';
my $noNT = 'bad NT length';
my $badNT = 'bad NT';
my $badGT = 'bad GT';
my $badDS1 = 'bad DS 1';
my $badDS2 = 'bad DS 2';


# just count the occurrences of the various strings that encode the genotypes
# e.g. 0/0:0:A|A

my $markers_to_analyze = shift // 100;
my $marker_count = 0;
my $accession_count = 0;
my %ok_pattern_count = ();
my %bad_pattern_count = ();
my @acc_badcounts = ();
my @marker_badcounts = ();
my @accids = ();
while (<>) {
  next if(/^\s*##/);
  if(/^\s*#/){
    my @acccols = split("\t", $_);
    $accession_count = scalar @acccols - 9;
    last;
  }
}
  while(<>){
  my @cols = split("\t", $_);
  my $ref = $cols[3];
  my $alt = $cols[4];
  #print "$ref $alt \n";
  my $refaltalleles = "$ref $alt";
  @cols = @cols[8..$#cols];
  my $format_str = shift @cols;
  my $n_colons_expected = () = $format_str =~ /:/g;
  for (my $i_acc=0; $i_acc < scalar @cols; $i_acc++) {
    my $s = $cols[$i_acc];
    # print "$s\n"; sleep(1);
    $s =~ s/\s+$//;
    my $p = "$refaltalleles   $s";
    my $p_bad = p_bad($p, $n_colons_expected);
    if($p_bad eq $OK){
      $ok_pattern_count{$p}++;
    }else{
      $bad_pattern_count{"$p  $p_bad"}++;
       $acc_badcounts[$i_acc]++;
      $marker_badcounts[$marker_count]++;
    }
  }
  $marker_count++;
  last if($marker_count >= $markers_to_analyze);
}
my $ok_count = sum(values %ok_pattern_count) // 0;
my $bad_count = sum(values %bad_pattern_count) // 0;
print "# total count: ", $ok_count + $bad_count, "  ok: $ok_count  bad: $bad_count \n";

@marker_badcounts = map($_ // 0, @marker_badcounts);
print "markers: $marker_count  markers with > 0 bad data: ", sum(map(($_ > 0)? 1 : 0, @marker_badcounts)), "\n";
@acc_badcounts = map($_ // 0, @acc_badcounts);
print "accessions: $accession_count  accessions with > 0 bad data: ", sum(map(($_ > 0)? 1 : 0, @acc_badcounts)), "\n";

print "Ok data:\n";
while (my ($p, $c) = each %ok_pattern_count) {
    print "$p   $c\n";
} print "\n";

print "Bad data:\n";
while (my ($p, $c) = each %bad_pattern_count) {
    print "$p   $c\n";
}

print "# total ok count,  bad count: $ok_count  $bad_count ";
printf ("(%5.3f %%)\n", 100*$bad_count/($ok_count + $bad_count));

sub p_bad{
  my $p = shift;
  my $n_colons = shift;
  my $RA = substr($p, 0, 3);
  my $GTDSNT = substr($p, 6);
  
  my $colon_count = () = $GTDSNT =~ /:/g;
  return $bad_colon_count if($colon_count != $n_colons);
  my ($GT, $DS, $NT) = split(":", $GTDSNT);
  return $noGT if((length $GT) ne 3);
  return $noDS if((length $DS) ne 1);
  return $noNT if((length $NT) ne 3);
  my @nts = ('A', 'C', 'G', 'T');
  for my $nt1 (@nts){
    for my $nt2 (@nts){
      next if($nt2 eq $nt1);
      if($RA eq "$nt1 $nt2"){
	my $NT1 = substr($GTDSNT, -3, 1);
	my $NT2 = substr($GTDSNT, -1, 1);
	return $badNT if(($NT1 ne $nt1 and $NT1 ne $nt2) or ($NT2 ne $nt1 and $NT2 ne $nt2)); # one of the NTs is neither ref nor alt allele
	my $GT1 = substr($GT, 0, 1);
	my $GT2 = substr($GT, 2, 1);
	my $GT1ex = ($NT1 eq $nt1)? 0 : 1;
	my $GT2ex = ($NT2 eq $nt1)? 0 : 1;
	return $badGT if($GT1 ne $GT1ex  or  $GT2 ne $GT2ex); # GT field inconsistent with NT, given what ref and alt alleles are.
	my $DSex = 0;
	$DSex++ if($NT1 eq $nt2);
	$DSex++ if($NT2 eq $nt2);
	return $badDS1 if(abs($DSex - $DS) == 1);
	return $badDS2 if($DSex ne $DS);
      }
    }
  }
  return $OK;
}
