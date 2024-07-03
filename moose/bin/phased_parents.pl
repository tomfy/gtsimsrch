#!/usr/bin/perl -w
use strict;


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
use Marker;

my $max_markers = 5000000;

my @acc_ids = ();
my %acc_chroms = (); # keys: acc ids, values: arrayref of chromosome object;
my %chrom_numbers = ();

while (my $line = <>) {
  next if($line =~ /^\s*##/);
  if ($line =~ /^\s*#/) {
    @acc_ids = split(" ", $line);
    @acc_ids = @acc_ids[9..$#acc_ids]; #
      # print "XXX ", scalar @acc_ids, "\n";
    last;
  }
}
# print 'YYY ',  scalar @acc_ids, "\n";
for my $acc_id (@acc_ids) {
  $acc_chroms{$acc_id} = [];
}

my $marker_count = 0;
while(my $line = <>) # process the first marker in the file
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
  while(my($i, $accid) = each @acc_ids){
    # print "i: $i  accid:  $accid\n";
    if(!exists $acc_chroms{$accid}->[$i_chrom]){ # 
      $acc_chroms{$accid}->[$i_chrom] = Chromosome->new({i_chrom => $i_chrom, genotypes => []});
    }
    my $the_chrom = $acc_chroms{$accid}->[$i_chrom];
    # print ref $the_chrom, "\n";
    # print $the_chrom->i_chrom(), "\n";
    # print scalar @{$the_chrom->genotypes()}, "\n";

    my $the_pgt = '4'; # 4 means missing data
    my $gt_info = $marker_gts[$i];
    my @fields = split(':', $gt_info);
    if($fields[0] eq '0|0'){
      $the_pgt = 0;
    }elsif($fields[0] eq '0|1'){
      $the_pgt = 1;
    }elsif($fields[0] eq '1|0'){
      $the_pgt = 2;
    }elsif($fields[0] eq '1|1'){
      $the_pgt = 3;
    }
    my $the_marker = Marker->new({# i_chrom => $i_chrom, id => $mrkr_id, position => $position, 
				  pgt => $the_pgt});
  $the_chrom->add_marker($the_marker);
  }
  $marker_count++;
  last if($marker_count >= $max_markers);
}

while(my($i, $accid) = each @acc_ids){
 # print "$i  $accid  ", scalar @{$acc_chroms{$accid}}, "\n";
  while (my ($i_chrom, $a_chrom) = each @{$acc_chroms{$accid}}){
    next if($i_chrom == 0);
  
 #   print "     ", $a_chrom->i_chrom(), "  ", join(' ', map($_->pgt, @{$a_chrom->genotypes()})), "\n";
  }
}
#sleep(2);

while(my ($i, $aid1) = each @acc_ids){
  my $pgts1 = $acc_chroms{$aid1}->[1]->genotypes();
  # print "XXX: ", join(" ", map($_->pgt(), @$pgts1)), "\n";
  # sleep(1);
#  while(my ($j, $aid2) = each @acc_ids){
  for (my $j = 0; $j < scalar @acc_ids; $j++){
    next if($j == $i);
    my $aid2 = $acc_ids[$j];
    my $pgts2 = $acc_chroms{$aid2}->[1]->genotypes();
    my ($Do, $So, $XA, $XB, $par1count, $par2count) = analyze_pgts_pair($pgts1, $pgts2);
    my $hgmr_denom = $Do + $So;
   print "$i $j  $aid1  $aid2   ", " $Do $hgmr_denom  ", ($hgmr_denom > 0)? $Do/$hgmr_denom : '-1', "   $XA $XB    $par1count $par2count ", $par1count + $par2count, "\n";
  }
}

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
  my $last_phaseA = undef;
  my $last_phaseB = undef;
  my $switch_countA = 0;
  my $switch_countB = 0;
  my $par_1count = 0;
  my $par_2count = 0;
 # print join('', map($_->pgt(), @$pgts1)), "\n";
 # print join('', map($_->pgt(), @$pgts2)), "\n";
  while (my($i, $pgt1) = each @$pgts1) {
    my $par_pgt = $pgt1->pgt();
    my $prog_pgt = $pgts2->[$i]->pgt();
    # print "$par_pgt  $prog_pgt \n";
    if ($par_pgt != 4  and  $prog_pgt != 4) {
      my $progA = ($prog_pgt == 0 or $prog_pgt == 1)? 0 : 1; # 0 <-> prog A has 0 (ref)
      my $progB = ($prog_pgt == 0 or $prog_pgt == 2)? 0 : 1; # 0 <-> prog B has 0 (ref)
      if ($par_pgt == 0) { # par: 0|0
	if ($prog_pgt == 0) { # prog: 0|0
	  $Shom++;
	} elsif ($prog_pgt == 3) { # prog: 1|1
	  $Dhom++;
	}
      } elsif ($par_pgt == 1) {  # par 0|1

	$par_1count++;
	if ($progA == 0) {
	  $phaseA = 0;
	  if (defined $last_phaseA  and  $phaseA != $last_phaseA) {
	    $switch_countA++;
	  }
	    $last_phaseA = $phaseA;
	} elsif ($progA == 1) {
	  $phaseA = 1;
	  if (defined $last_phaseA  and  $phaseA != $last_phaseA) {
	    $switch_countA++;
	  }
	    $last_phaseA = $phaseA;
	}
	if ($progB == 0) {
	  $phaseB = 0;
	  if (defined $last_phaseB  and  $phaseB != $last_phaseB) {
	    $switch_countB++;
	  }
	    $last_phaseB = $phaseB;
	} elsif ($progB == 1) {
	  $phaseB = 1;
	  if (defined $last_phaseB  and  $phaseB != $last_phaseB) {
	    $switch_countB++;
	  }
	    $last_phaseB = $phaseB;
	}
#			print "X $par_pgt  $prog_pgt   $progA $progB  $last_phaseA $last_phaseB   $switch_countA $switch_countB \n";

      } elsif ($par_pgt == 2) {  # par 1|0
	$par_2count++;
	if ($progA == 0) {
	  $phaseA = 1;
	  if (defined $last_phaseA  and  $phaseA != $last_phaseA) {
	    $switch_countA++;
	    $last_phaseA = $phaseA;
	  }
	} elsif ($progA == 1) {
	  $phaseA = 0;
	  if (defined $last_phaseA  and  $phaseA != $last_phaseA) {
	    $switch_countA++;
	    $last_phaseA = $phaseA;
	  }
	}
	if ($progB == 0) {
	  $phaseB = 1;
	  if (defined $last_phaseB  and  $phaseB != $last_phaseB) {
	    $switch_countB++;
	  }
	    $last_phaseB = $phaseB;
	} elsif ($progB == 1) {
	  $phaseB = 0;
	  if (defined $last_phaseB  and  $phaseB != $last_phaseB) {
	    $switch_countB++;
	  }
	    $last_phaseB = $phaseB;
	}
      }elsif( $par_pgt == 3) { # par: 1|1
	if ($prog_pgt == 0) { # prog: 0|0
	  $Dhom++;
	} elsif ($prog_pgt == 3) { # prog: 1|1
	  $Shom++;
	}
      }
      #  print STDERR "$par_pgt  $progA $progB  " , $phaseA // '-', "  $switch_countA  ", $phaseB // '-', "  $switch_countB\n";
    }
  }
  return ($Dhom, $Shom, $switch_countA, $switch_countB, $par_1count, $par_2count);
}
