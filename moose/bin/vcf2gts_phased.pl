#!/usr/bin/perl -w
use strict;
use List::Util qw(min max sum shuffle);

# read in a vcf file
# if it has phased genotypes, analyze to find
# likely parent-progeny pairs.

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

# analyze a vcf file and use the phase information to
# test pedigrees.
# genotypes fields of the vcf file must have
# phased information as first subfield
# if 3rd subfield has probabilities for the 3 possible dosages
#    the max of these 3 will be required to be > $min_gt_prob
#    otherwise the genotypes is considered to be missing.
# for example:  1|0:1:0,1,0
# In this case the probability of the genotype is 1, so it is retained.

# usage: phased_parents.pl -ped <pedigree filename> -vcf < <vcf filename> [-out <output filename>] [-rand <n random parents>] [-min_prob <min gt probability>] [-rev (to do in rev direction also>]



#my $pedigree_file = undef;
my $vcf_file = undef;
my $output_file = undef;
#my $rand_parents_to_do = 0; # for each accession with pedigree, also choose this many other pairs of accessions at random to test as parents.
my $min_gt_prob = 0.9;	    # for 
#my $do_reverse = 0;
my $max_markers = 5000000; # to limit to a smaller number of markers for speed in testing
#my $use_pedigrees = 0;
my $help = 0;
#my $remove_length_one = 0;

GetOptions(
	 #  'pedigree_file=s' => \$pedigree_file,
	   'vcf_file=s' => \$vcf_file,
	   'output_file=s' => \$output_file,
	   #'rand_parents|n_rand_parents=i' => \$rand_parents_to_do,
	   'min_prob|gt_min_prob=f' => \$min_gt_prob,
	   #'reverse!' => \$do_reverse,
	   'marker_limit=i' => \$max_markers,
	   'help!' => \$help,
	   #'remove_length_one!' => \$remove_length_one,
	  );

if(!defined $vcf_file  or  $help){
  print_help();
  exit;
}

if(!defined $output_file){
  $output_file = $vcf_file . '.pdsgs';
}

open my $fhout, ">", "$output_file";
my $run_parameter_info = 
  "# vcf file: $vcf_file\n".
  "# output file: $output_file\n".
  "# min gt probability: $min_gt_prob\n";
print STDOUT $run_parameter_info;
print $fhout $run_parameter_info;

##############################################
# read and store vcf file ####################
##############################################

my @acc_ids = (); # acc ids in vcf file
my %acc_chrom_objs = (); # keys: acc ids, values: arrayref of chromosome object;
my %chromnumber_x = (); # x for exists. keys are chrom number which occur in vcf, values: 1
my @marker_ids = ();
my %mrkrid_chrom = ();
my @chrom_numbs = (); # in the order they occur in vcf file.
my %chrom_markerid = (); # keys: chrom numbers, values: array ref of marker ids on chrom (in order as in vcf file)

open my $fhvcf, "<", "$vcf_file";
while (my $line = <$fhvcf>) {	# read accession ids from vcf file
  next if($line =~ /^\s*##/);
  next if($line =~ /^\s*$/);
  if ($line =~ /^\s*#/) {
    @acc_ids = split(" ", $line);
    @acc_ids = @acc_ids[9..$#acc_ids]; #
    print STDOUT "# number of accession ids in vcf file: ", scalar @acc_ids, "\n";
    last;
  }
}

for my $acc_id (@acc_ids) {
  $acc_chrom_objs{$acc_id} = [];
}

# read the vcf file and store the phased genotype information
my $marker_count = 0;
while (my $line = <$fhvcf>) # each pass through this loop processes one marker, all accessions
  {
    my @marker_gts = split(" ", $line); # these are the (phased) genotypes for the various accessions, this marker
    # things like  0|0:0:1,0,0
    # i.e. colon separated fields. here we require the first of these to be '0|0' '0|1' '1|0', '1|1'
    my $i_chrom = shift @marker_gts;
    $chromnumber_x{$i_chrom} = 1;
    push @chrom_numbs, $i_chrom;
    my $position = shift @marker_gts;
    my $mrkr_id = shift @marker_gts;
    my $ref = shift @marker_gts;
    my $alt = shift @marker_gts;
    my $qual = shift @marker_gts;
    my $filter = shift @marker_gts;
    my $info = shift @marker_gts;
    my $format = shift @marker_gts;
    push @marker_ids, $mrkr_id;
    $mrkrid_chrom{$mrkr_id} = $i_chrom;
      if (!exists $chrom_markerid{$i_chrom}) {
	$chrom_markerid{$i_chrom} = [$mrkr_id];
      } else {
	push @{$chrom_markerid{$i_chrom}}, $mrkr_id;
      }
    while (my($i, $accid) = each @acc_ids) { # loop over the accessions
      if (!exists $acc_chrom_objs{$accid}->[$i_chrom]) { # check whether Chromosome obj. exists for this accession and chromosome number
	$acc_chrom_objs{$accid}->[$i_chrom] = Chromosome->new({i_chrom => $i_chrom, genotypes => []});
      }
    
      my $the_chrom = $acc_chrom_objs{$accid}->[$i_chrom];
      my $the_pgt = 'X';	# X means missing data
      my $gt_info = $marker_gts[$i];
      my @fields = split(':', $gt_info);
      my @gps = split(',', $fields[2]); # est probabilities for different 
      if (max(@gps) >= $min_gt_prob) { # good quality data for this accession/marker combo
	if ($fields[0] eq '0|0') {
	  $the_pgt = 0;
	} elsif ($fields[0] eq '0|1') {
	  $the_pgt = 1;
	} elsif ($fields[0] eq '1|0') {
	  $the_pgt = -1; # 2;
	} elsif ($fields[0] eq '1|1') {
	  $the_pgt = 2; # 3;
	}
      }
      $the_chrom->add_genotype($the_pgt);
    }
    $marker_count++;
    printf STDOUT "# Markers read from vcf file so far: $marker_count \n" if($marker_count % 200 == 0);
    last if($marker_count >= $max_markers);
  }
printf STDOUT "# Done storing markers. Markers stored:  $marker_count \n";
close $fhvcf;
###########################################
# done storing info in vcf file. ##########
###########################################

# my $c = getc();

my @sorted_chrom_numbers = sort {$a <=> $b} keys %chromnumber_x; # 
print STDERR "Chromosome numbers: [", join(", ", @sorted_chrom_numbers), "]\n";
#sleep(4);
#my @sorted_chroms = 

my $marker_id_string = "MARKER ";
my $chromosome_number_string = "CHROMOSOME ";
#for my $mrkrid (@marker_ids){ # this is in the order they occur in vcf file.
#  print $fhout $mrkrid_chrom{$mrkrid}, " ";
#} print $fhout "\n";

my $first_accid = $acc_ids[0];
my $the_chroms = $acc_chrom_objs{$first_accid};
while( my ($i, $achrom) = each @$the_chroms){
  next if($i == 0  or  (!ref($achrom)));

  my $the_markerids = $chrom_markerid{$achrom->i_chrom()};
    print STDERR "XTZ i  ", $achrom->i_chrom(), "  ", scalar @$the_markerids, "\n";
  $marker_id_string .= " " . join(" ", @$the_markerids);
  my @cs = map($achrom->i_chrom(), @$the_markerids);
  $chromosome_number_string .= " " . join(" ", @cs);
}
print $fhout $marker_id_string, "\n";
print $fhout $chromosome_number_string, "\n";

#sleep(4);
my $chrom_line_done = 0; # set to 1.
for my $accid (@acc_ids){
  # 0 -> 0,  1 -> 1,  2 -> -1,  3 -> 2
  my $the_chroms = $acc_chrom_objs{$accid}; # $the_chroms is an array ref of Chromosome objs.
  print STDERR "$accid   number of chroms: ", scalar @$the_chroms, "\n";
  printf $fhout "%24s  ", $accid;
  while( my ($i, $achrom) = each @$the_chroms){
    print STDERR "$accid  $i  [", ref($achrom), "]\n";
    next if($i == 0  or  (!ref($achrom)));
    die "chrom numbers don't match.\n" if($i != $achrom->i_chrom());
  
    #sleep(1);
    my $the_gts = $achrom->genotypes();
    # printf $fhout " ", join(" ", @$the_gts);
    for my $agt (@$the_gts){
      if($agt eq 'X'){
	printf $fhout "   X";
      }else{
	printf $fhout " %3d", $agt;
      }
    }
  }
  print $fhout "\n";
}
exit;

sub print_help{ 
    print "Help\n";
}

