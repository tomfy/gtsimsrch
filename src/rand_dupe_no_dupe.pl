#!/usr/bin/perl -w
use strict;
use List::Util qw(shuffle);

my $genofile = shift;		# .ped format
my $phenofile = shift;         # id id  phenotype
my $alleles_file = shift; # typical line: S3_143567 A T - needed for adding errors to duplicates.

# my $serial_number = shift // 0;
my $ufilename = shift;
my $udfilename = shift;

my $nunique = shift // 300;
my $nduplicate = shift // 150;
#print STDERR "$nunique  $nduplicate \n";
#sleep(1);
my $pnoise = shift // 0.01; # shift // 0.1;

#%accid_ped_line = ();

my $u_noisy = 1;

# store genotypes from file
my %accid_gts = ();
open my $fhin, "<", "$genofile";
while (<$fhin>) {
  if (/^\s*(\S+)(.*)$/) {
    my ($accid, $gtstr) = ($1, $2);
    $accid_gts{$accid} = $gtstr;
  }
}
close $fhin;

# store phenotypes from file
my %accid_blup = ();
open $fhin, "<", "$phenofile";	#
while (<$fhin>) {
  my ($ida, $idb, $blup) = split(" ", $_);
  #print STDERR "$ida $idb $blup\n";
  $accid_blup{$ida} = $blup if(exists $accid_gts{$ida}); # don't store if no genotypes
}
close $fhin;
print "# number of accession with both genotypes and phenotypes: ", scalar keys %accid_blup, "\n";
#exit;

# store alleles from file
my @allelepairs = ();
open $fhin, "<", "$alleles_file";
while (<$fhin>) {
  my ($id, $A1, $A2) = split(" ", $_);
  push @allelepairs, ($A1, $A2);
}
close $fhin;

# if random subset requested, shuffle and take first $nunique
my @ids;
if($nunique < scalar keys %accid_gts){
  @ids = shuffle(keys %accid_gts);
  @ids = @ids[0..$nunique-1];
}else{
  @ids = keys %accid_gts;
}

# get unique ped file strings (with noise added)
# and corresponding phenotype (blup) strings
my @orig_u_ped_strings = ();
my @orig_u_blup_strings = ();

my @u_ped_strings = ();
my @u_blup_strings = ();

my $unique_gtstr = '';
my $unique_blupstr = '';
  for my $id (@ids){
  my $ped_line =  "$id  " . $accid_gts{$id} . "\n";
  push @orig_u_ped_strings, $ped_line;

  my @cols = split(" ", $ped_line);
  my $id = $cols[0];
  #print STDERR "#last col index:  ", $#cols, "\n";
  my @noisy_cols = ($u_noisy)? @{ add_errors_to_ped_line(\@cols, \@allelepairs, $pnoise) } :
    @cols[6..$#cols];
  my $noisy_u_ped_string = "$id $id  0 0 0 -9 " . join(" ", @noisy_cols) . "\n";
  
#  print STDERR "Orig. AGMRs:  ", agmr($ped_line, $noisy_u_ped_string), "\n";
  # print $ped_line;
  # print $noisy_u_ped_string;
  # exit;
  push @u_ped_strings, $noisy_u_ped_string;

  push @orig_u_blup_strings,  "$id $id  " . $accid_blup{$id} . "\n";
  # add noise to blups here - not implemented 
  push @u_blup_strings,  "$id $id  " . $accid_blup{$id} . "\n";
}

#print STDERR scalar @u_ped_strings, "\n"; sleep(1);

# print the subset with noise added
#my $ufilename = "u" . $nunique . "_" . $serial_number;
my $ugtfilename = $ufilename . ".ped";
my $uphfilename = $ufilename . ".phen";
open my $fhout, ">", $ugtfilename;
print $fhout join('', @u_ped_strings);
close $fhout;
open $fhout, ">", $uphfilename;
print $fhout join('', @u_blup_strings);
close $fhout;

#@ids = @ids[0..$nunique-1];

# @ids is the array of a subset of unique acc ids
my @dupe_ped_lines = ();
my @dupe_blup_lines = ();
my $dupe_number = 1;
for my $idupe (1..$nduplicate) {
  #print STDERR "Working on duplicate $dupe_number.\n";
  my $irand = int(rand($nunique));
  #print STDERR "$irand  ", scalar @u_ped_strings, "\n";
  my $ped_line = $orig_u_ped_strings[$irand];
  my $blup_str = $u_blup_strings[$irand];
#  print STDERR "$irand  ", $blup_str // 'undef ', "\n";
  my ($x, $y, $blup) = split(" ", $blup_str);
  my @cols = split(" ", $ped_line);
  my $id = $cols[0] . '_D' . $dupe_number;
  #print STDERR "#last col index:  ", $#cols, "\n";
  my @noisy_cols = @{ add_errors_to_ped_line(\@cols, \@allelepairs, $pnoise) };
  #$noisy_cols[0] = $id;
  #$noisy_cols[1] = $id;
  my $dupe_ped_line = "$id $id  0 0 0 -9 " . join(" ", @noisy_cols) . "\n";
 # print STDERR "Dupe. AGMRs:  ", agmr($ped_line, $dupe_ped_line), "\n";
  push @dupe_ped_lines, $dupe_ped_line;
  push @dupe_blup_lines, "$id $id  $blup\n";
  $dupe_number++;
}

#my $dgtfilename = "orig_plus_dupes.ped";
#my $dphfilename = "orig_plus_dupes.pheno";
#my $dfilename = "d" . $nduplicate . "_" . $serial_number;
my $dgtfilename = "d.ped";
my $dphfilename = "d.phen";
open $fhout, ">", $dgtfilename;
#print $fhout join('', @u_ped_strings);
print $fhout join('', @dupe_ped_lines);
close $fhout;
open $fhout, ">", $dphfilename;
#print $fhout join('', @u_blup_strings);
print $fhout join('', @dupe_blup_lines);
close $fhout;

#my $udfilename = "u" . $nunique . "d" . $nduplicate . "_" . $serial_number;
my $udgtfilename = $udfilename . ".ped";
my $udphfilename = $udfilename . ".phen";

system "cat $ugtfilename $dgtfilename  > $udgtfilename";
system "cat $uphfilename $dphfilename  > $udphfilename";

unlink $dgtfilename;
unlink $dphfilename;

sub add_errors_to_ped_line{	# switch some of the genotypes
  my $car = shift;		# array ref
  my $allelepairs = shift;
  my $pnoise = shift;
  #my @cols = split(" ", $ped_line);
  my @cols = @$car;
  @cols = @cols[6..$#cols];
  # my $id = $cols->[0];
  # my @cols = @cols[6..$#cols];
  #print STDERR "last col index: ", $#cols, "\n";
  for (my $i = 0; $i <= $#cols; $i += 2) {
    my $randnum = rand(1);
    #  print STDERR "$i  $randnum  $pnoise\n";
    if ($randnum < $pnoise) {	# replace with a random gt
      #    print STDERR "XXXXX\n";
      my $Ai = $cols[$i];
      my $Aj = $cols[$i+1];
      my $x = "$Ai $Aj ";
      my $A1 = $allelepairs->[$i];
      my $A2 = $allelepairs->[$i+1];
      #print STDERR "$i {$A1}{$A2}  {$Ai}\n" if(!defined $A1);
      die "A1 A2 Ai: $A1 $A2  $Ai\n" if($Ai ne $A1  and  $Ai ne $A2);
      $cols[$i] = ($Ai eq $A1)? $A2 : $A1;
      #$cols[$i+1] = (rand(1) < 0.5)? $A1 : $A2;
      # print STDERR "$x  ", $cols[$i], " ", $cols[$i+1], "\n";
    }
  }
  #return "$id  $id  0 0 0 -9  " . join(" ", @cols) . "\n";
  return \@cols;
}

sub agmr{
  my $p1 = shift;
  my $p2 = shift;
  my $agmr = 0;
  my @cols1 = split(" ", $p1);
  @cols1 = @cols1[6..$#cols1];
  my @cols2 = split(" ", $p2);
  @cols2 = @cols2[6..$#cols2];
  # print STDERR scalar @cols1, " ", scalar @cols2, "\n";
  for(my $i=0; $i<= $#cols1; $i+= 2){
    my $gt1 = $cols1[$i] . '_' . $cols1[$i+1];
    my $gt2 = $cols2[$i] . '_' . $cols2[$i+1];
    $agmr++ if($gt1 ne $gt2);
  }
  return $agmr;
}
  
