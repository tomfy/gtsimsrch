#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw(min max sum shuffle);
use File::Spec qw(splitpath);
use File::Basename 'dirname';

use Cwd 'abs_path';
my $bindir;
BEGIN {     # this has to go in Begin block so happens at compile time
  $bindir =
    dirname( abs_path(__FILE__) ) ; # the directory containing this script
}

# reads in a dosage_matrix file (e.g. output of vcf_to_gts)
# and a cluster file (output of clusterer)
# and optionally a pedigree file 
# outputs a dosage matrix file with duplicates removed, i.e.
# with the same format as duplicate_search input, but now with just one
# line representing each cluster of duplicates.
# also output a uniquified pedigree file if specify a pedigree file to be input.
# optionally ( -vote ) have the members of each group of duplicates vote on each genotype.

my $dosages_input_filename = undef; # must be specified 
my $pedigrees_input_filename = undef; # if specified, pedigree file to be uniquified
my $dosages_output_filename = undef; # if not specified prepend 'u_' to $dosages_input_filename
my $pedigrees_output_filename = undef; # if not specified prepend 'u_' to $pedigree_input_filename
my $cluster_filename = undef; # output of clusterer, with info on clusters

my $cluster_fraction = 0.0; # fraction of other cluster members to keep (aside from one representative which is always kept)
my $vote = 0;
my $missing_str = "X";
my $progeny_col = -3; # default is pedigree ids are in last 3 cols (progeny, Fparent, Mparent)
my $Fparent_col = -2;
my $Mparent_col = -1;
my $nskip = 0; # number of lines to skip at top of pedigree file.

GetOptions(
	   'dosages_input_file=s' => \$dosages_input_filename,
	   'dosages_output_file=s' => \$dosages_output_filename, # uniquified
	   'pedigrees_input_file=s' => \$pedigrees_input_filename,
	   'pedigrees_output_file=s' => \$pedigrees_output_filename, # uniquified
	   'cluster_filename=s' => \$cluster_filename,

# less useful options:
	   'fraction=f' => \$cluster_fraction,
	   'vote!' => \$vote,
	   'missing_str=s' => \$missing_str,
	  );
if(!defined $dosages_input_filename  or  !defined $cluster_filename){
if (!defined $dosages_input_filename) {
  print "Must specify dosages file name.\n";
}
if(!defined $cluster_filename){
  print "Must specify cluster filename (output of clusterer).\n";
}
  usage_message();
  exit;
}
my ($v, $dir, $dosages_filename_in) = File::Spec->splitpath( $dosages_input_filename );
if (!defined $dosages_output_filename) {
  $dosages_output_filename = "u_" . $dosages_filename_in;
  print STDERR "uniquified dosages file:    $dosages_output_filename \n";
}
###############################################################################
####    read in cluster information         ###################################
###############################################################################
open my $fh_clusters, "<", "$cluster_filename";
my %clusterid_repid = (); # keys: all ids of accessions in clusters; values: ids of corresponding representative accessions.
my %repid_clusterids = (); # keys: ids of cluster representatives; values: array refs with corresponding cluster ids.

while (my $line = <$fh_clusters>) { # each line is one cluster
  next if($line =~ /^\s*#/);
  my @cols = split(" ", $line);
  next if(scalar @cols  < 11); # cluster line should have 9 quality params plus at least 2 cluster ids.
  my $cluster_size = shift @cols;
  my $min_d = shift @cols;
  my $mean_d = shift @cols;
  my $max_d = shift @cols;
  my $min_intraextra_d = shift @cols;
  my $n_near1 = shift @cols;
  my $n_near2 = shift @cols;
  my $n_big_intra_d = shift @cols;
  my $min_cluster_member_degree = shift @cols; # min degree of cluster members

  my @cluster_accids = ();
  if (scalar @cols == $cluster_size) { # short format: just cluster ids
    @cluster_accids = @cols;
  } elsif (scalar @cols == (3*$cluster_size)) { # long format: each id followed by degree, and ___
    for (my $i = 0; $i < scalar @cols; $i += 3) {
      push @cluster_accids, $cols[$i];
    }
  } else {
    die "cluster size inconsistent. $cluster_size ", scalar @cols;
  }
  die "cluster size inconsistent. $cluster_size ", scalar @cluster_accids if($cluster_size != scalar @cluster_accids);

  my $rep_id = $cluster_accids[0];     # id of the representative of the cluster
  # print STDERR "repid $rep_id \n"; # sleep(1);
  $repid_clusterids{$rep_id} = \@cluster_accids;
  for my $a_cluster_id (@cluster_accids) {
    $clusterid_repid{$a_cluster_id} = $rep_id; # all accessions in clusters get stored
  }
}
close $fh_clusters;

###############################################################################
####   read in pedigree file, if specified.    ################################
###############################################################################
if(defined $pedigrees_input_filename){
  open my $fhpedin, "<", "$pedigrees_input_filename";
  if(!defined $pedigrees_output_filename){
    my ($v, $dir, $pedigrees_filename_in) = File::Spec->splitpath( $pedigrees_input_filename );
    $pedigrees_output_filename = "u_" . $pedigrees_filename_in;
    print STDERR "uniquified pedigrees file:  $pedigrees_output_filename \n";
  }
  # want col 1 to be leftmost, -1 rightmost
$progeny_col-- if($progeny_col > 0);
$Fparent_col-- if($Fparent_col > 0);
$Mparent_col-- if($Mparent_col > 0);
open my $fhpedout, ">", "$pedigrees_output_filename";
for(1..$nskip){ <$fhpedin>; } # skip first nskip lines.

my $pedigrees_in_count = 0;
my $progeny_NA_count = 0;
my $duplicate_progeny_count = 0;
my $pedigrees_out_count = 0;

my %uniqprogeny = ();
while (my $line = <$fhpedin>) {
  next if($line =~ /^\s*#/);
  $pedigrees_in_count++;
  $line =~ s/\s+$//;	# remove newline, any other whitespace at end.
  my @cols = split(" ", $line);
  my $progeny_id = $cols[$progeny_col];

  if ($progeny_id ne 'NA') { # output, replacing ids by ids of the cluster representatives if in a cluster
    if (exists $uniqprogeny{$progeny_id}) { # the progeny is a duplicate of a previous one.
      $duplicate_progeny_count++;
    } else {
      $uniqprogeny{$progeny_id}++;
      $pedigrees_out_count++;
      $progeny_id = $clusterid_repid{$progeny_id} // $progeny_id;
      my $Fpar_id = $clusterid_repid{$cols[$Fparent_col]} // $cols[$Fparent_col]; # replace id with rep id if F parent is a duplicate group member
      my $Mpar_id = $clusterid_repid{$cols[$Mparent_col]} // $cols[$Mparent_col]; # replace id with rep id if M parent is a duplicate group member
      print $fhpedout "$progeny_id  $Fpar_id  $Mpar_id\n";
    }
  } else {
    $progeny_NA_count++;
  }
}
close $fhpedin;
} # end if defined $pedigrees_input_filename

################################################################################
####   read in dosages  ########################################################
####   store id and dosages of cluster members    ##############################
####   output id and dosages of others       ###################################
################################################################################
print STDERR "before storing dosages\n";

# store individual lines of dosage file in hash
open my $fh_dosages, "<", "$dosages_input_filename";
open my $fhout, ">", "$dosages_output_filename" or die "Couldn't open $dosages_output_filename for writing.\n";

# store ids and genotypes of clusters (size >= 2), and output singletons.
my %id_gts  = ();  # key: ids; value: array ref of dosages (0,1,2,NA) 
my $first_line = <$fh_dosages>;
print $fhout $first_line;
while (my $line = <$fh_dosages>) {
  if ($line =~ /^\s*(#|MARKER|CHROMOSOME)/) {
    print $fhout $line;
    next;
  }
  my @cols = split(" ", $line);
  my $id = shift @cols;
  if (exists $clusterid_repid{$id}) { # this accession belongs to a cluster
    $id_gts{$id} = \@cols; # value is array ref with gts as in dosage file 
  } else { # accession not in a cluster, just output the line just read in.
    print $fhout $line;
  }
}
close($fh_dosages);
print STDERR "after storing dosages\n";
###############################################################################

my @duplicate_lines = ();

###############################################################################
####   output uniquified dosage matrix              ###########################
###############################################################################
if ($vote  and  $cluster_fraction == 0) { # cluster members vote, and then output just representative id with 'elected' dosages.
  print STDERR "duplicate group genotypes determined by voting.\n";
  while (my($representative_id, $cluster_accids) = each %repid_clusterids) {
    my $elected_gts = vote($cluster_accids, \%id_gts);
    print $fhout "$representative_id  ", join(" ", @$elected_gts), "\n";
  }
} else {
  print STDERR "duplicate group representative accessions are those  with least missing data.\n";
  while (my($representative_id, $cluster_accids) = each %repid_clusterids) {
    print $fhout "$representative_id  ", join(" ", @{$id_gts{$representative_id}}), "\n"; # output representative and its dosages.
  }
}

###############################################################################
# if $cluster_fraction > 0, output specified fraction of duplicates:
@duplicate_lines = shuffle @duplicate_lines;
my $n_duplicates_to_output = int($cluster_fraction * scalar @duplicate_lines + 0.5);
for my $i (1..$n_duplicates_to_output) {
  print $fhout $duplicate_lines[$i];
}
close $fhout;

##################################################################################
###   end of main   ##############################################################
##################################################################################


sub vote{ # for clusters with several members, determine 
  my $cluster_ids = shift; # array ref holding the ids of the accessions in the cluster.
  my $id_dosages = shift; # hash ref. key: id, value: array ref of dosages for this id.
  my $first_id = $cluster_ids->[0];
  my $first_dosages = $id_dosages->{$first_id};
  my $n_markers = scalar @$first_dosages;
  my @elected_dosages = (0) x $n_markers;
  my $cluster_size = scalar @$cluster_ids;

  if ($cluster_size > 2) {
    my @marker_votes = ();;
    for (1..$n_markers) {
      push @marker_votes, [0, 0, 0, 0];
    }
    for my $an_id (@$cluster_ids) { # loop over accessions in cluster
      my $dosages = $id_dosages->{$an_id}; # array ref of dosages for $an_id
      while (my($i, $d) = each @$dosages) {
	if ($d eq $missing_str) {
	  $d = 3;
	}
	$marker_votes[$i]->[$d]++;
      }
    }
    while (my($i, $mv) = each @marker_votes) {
      my $e = $missing_str;
      for my $j (0..2) {
	my $votes = $mv->[$j]; # the number of votes for dosage = $j for marker $i
	if (2*$votes > $cluster_size) { # must have > 50% for one dosage, or $missing_str
	  $e = $j;
	  last;
	}
      }
      $elected_dosages[$i] = $e;
    }
     my $X_count = sum(map( (($_ eq 'X')? 1 : 0), @elected_dosages)); #  =~ tr/X/X/;
  } else {			# cluster size == 2
    my $second_id = $cluster_ids->[1];
    while (my($i, $d1) = each @$first_dosages) {

      my $d2 = $id_dosages->{$second_id}->[$i];
      if ($d1 eq $d2) {
	$elected_dosages[$i] = $d1;
      } elsif ($d1 eq $missing_str) { # d1 is missing, use d2 (which is not missing)
	$elected_dosages[$i] = $d2;
      } elsif ($d2 eq $missing_str) {	# d2 is missing, use d1 (which is not missing)
	$elected_dosages[$i] = $d1;
      } else {		    # neither is missing, and they disagree; use $missing_str
	$elected_dosages[$i] = $missing_str;
      }
    }
  } # end cluster size == 2 branch
  return \@elected_dosages;
}

sub usage_message{
  print "Basic usage: uniquify -dosages_in <input filename> -cluster <cluster filename> [-pedigrees_in <pedigree file>].\n";
  print "Options:\n";
  print "  -dosages_input_file  <dosage file to uniquify. Required>\n";
  print "  -dosages_output_file <name for  uniquified dosage file. default: prepend 'u_' to input filename.\n";
  print "  -pedigrees_input_file <pedigree file to uniquify. optional>\n";
  print "  -pedigrees_output_file <name for uniquified pedigree file.\n";
  print "  -cluster_filename <input file specifying the clusters. Required.\n";
}

# ************************************************
# sub least_missing_data_id{ # not used, now this is done in clusterer
#   my $cluster_ids = shift;	# array ref
#   my $id_genotypes = shift;	# hash ref, values: array ref of gts;
#   my $least_md_count = 1000000000000;
#   my $least_md_id = undef;
#   my $first_md_count = undef;
#   my $marker_count = undef;
#   while(my($i, $id) = each @$cluster_ids) {
#     if (exists $id_genotypes->{$id}) {
#       my @gts = @{$id_genotypes->{$id}};
#       my $md_count = sum(map( (($_ eq 'X')? 1 : 0), @gts));
#       $first_md_count = $md_count if($i == 0);
#       $marker_count = scalar @gts;
#       if($md_count < $least_md_count){
# 	$least_md_count = $md_count;
# 	$least_md_id = $id;
#       }
#     }else{
#       print STDERR "# Genotypes not available for $id \n";
#     }
#   }
#   print STDERR "AAAAA  $marker_count   $least_md_id  $least_md_count ", $cluster_ids->[0], " $first_md_count \n";
#   return $least_md_id;
# }

