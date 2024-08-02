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

# runs duplicatesearch, then clusterer, and outputs a file
# with the same format as duplicatesearch input, but now with just one
# line representing each cluster of duplicates.
# optionally ( -vote ) have the members of each group of duplicates vote on each genotype.
# optionally (  -filter_out <afilename>  ) output a dosage matrix file which is filtered as well as uniquified.

my $input_dosages_filename = undef;
#my $do_remove_bad_accessions = shift // m
my $max_acc_missing_data_fraction = 0.5;
my $max_marker_missing_data_fraction = 0.25;
my $min_maf = 0.08;
my $output_dosages_filename = undef;
my $max_distance = undef;
my $cluster_max_distance = 'auto';
my $cluster_fraction = 0.0; # fraction of other cluster members to keep (aside from one representative which is always kept)
my $vote = 0;
my $missing_str = "X";
my $filter_out = 1; # output the filtered set of markers, as well as removing duplicates
my $rng_seed = undef; # default is to let duplicatesearch get seed from clock
my $phenotype_filename = undef;
my %accid_phenotype = ();

GetOptions(
	   'input_file=s' => \$input_dosages_filename,
	   'output_file=s' => \$output_dosages_filename,
	   'acc_max_md_fraction=f' => \$max_acc_missing_data_fraction,
	   'marker_max_md_fraction=f' => \$max_marker_missing_data_fraction,
	   'distance_max=f' => \$max_distance,
	   'min_maf|maf_min=f' => \$min_maf,
	   'cluster_max_distance=f' => \$cluster_max_distance,
	   'fraction=f' => \$cluster_fraction,
	   'vote!' => \$vote,
	   'missing_str=s' => \$missing_str,
	   'filterout|filter_out!' => \$filter_out,
	   'seed|rng_seed=i' => \$rng_seed,
	   'phenotype_filename=s' => \$phenotype_filename, # id id phenotype ; col 3 is a quantitative phenotype, cols 1 and 2 are the same accession id twice
	  );

if (!defined $input_dosages_filename) {
  print "Must supply input file name.\n";
  usage_message();
  exit;
}
my ($v, $dir, $dosages_filename) = File::Spec->splitpath( $input_dosages_filename );
if (!defined $output_dosages_filename) {
  $output_dosages_filename = $dosages_filename . "_filtered" if($filter_out);
  $output_dosages_filename .= "_duplicates_removed";
  print STDERR "$dir    $output_dosages_filename \n";
}
#exit;

#  if a phenotype file is specified, store its id-phenotype pairs.
if(defined $phenotype_filename){
  open my $fhphen, "<", "$phenotype_filename";
  while( <$fhphen> ){
    next if(/^\s*#/);
    my ($id1, $id2, $pheno) = split(" ", $_); # id1 and id2 should be the same
    $accid_phenotype{$id1} = $pheno;
  }
  close $fhphen;
}


# open my $fhin, "<", "$input_dosages_filename" or die "couldn't open $input_dosages_filename for reading.\n";
my $cleaned_dosages_filename = $input_dosages_filename . "_cleaned";

####   remove accessions with excessive missing data   #########################

# if (0) {
#   open my $fhout, ">", "$cleaned_dosages_filename";
#   my $n_markers = undef;
#   my $n_bad_accessions = 0;
#   print STDERR "# Removing accessions with too much missing data.\n";
#   while (my $line_in = <$fhin>) {
#     if ($line_in =~ /^\s*#/) {
#       # print STDERR "print comment line.\n";
#       print $fhout $line_in;
#     } elsif ($line_in =~ /^MARKER/) {
#       #  print STDERR "print marker ids line.\n";
#       print $fhout $line_in;
#       my @markers = split(" ", $line_in);
#       $n_markers = scalar @markers  - 1;
#     } else {
#       if($line_in =~ /^\s*(\S+)/){
# 	my $accession_id = $1;
#       die "File lacks line with MARKER and marker ids\n" if(! defined $n_markers);
#       my $n_bad = () = $line_in =~ /\sNA/g;
#       if ($n_bad/$n_markers <= $max_acc_missing_data_fraction) {
# 	print $fhout $line_in;
#       } else {
# 	$n_bad_accessions++;
# 	print STDERR "Removing accession $accession_id, which has missing data for $n_bad markers.\n";
#       }
# }
#     }
#   }
#   close $fhout;
#   print STDERR "# $n_bad_accessions accessions eliminated due to excessive missing data.\n";
#   print "# $n_bad_accessions accessions eliminated due to excessive missing data (>" ,
#     int($max_acc_missing_data_fraction*100 + 0.5), "\%)\n";
# } else {
# system "~/gtsimsrch/src/bad_accessions_begone.pl -i $input_dosages_filename -o $cleaned_dosages_filename -m $max_acc_missing_data_fraction";
# }

print STDERR "dosages file with high-missing data accessions removed: $cleaned_dosages_filename \n";
#exit;

my $duplicatesearch_command = "duplicatesearch  -in $input_dosages_filename ";
$duplicatesearch_command .= " -distance $max_distance " if(defined $max_distance);
$duplicatesearch_command .= " -accession_max_missing_data  $max_acc_missing_data_fraction ";
$duplicatesearch_command .= " -maf_min $min_maf ";
$duplicatesearch_command .= "-marker_max_missing_data $max_marker_missing_data_fraction " if(defined $max_marker_missing_data_fraction);
$duplicatesearch_command .= " -seed $rng_seed " if(defined $rng_seed);
my $filtered_dosages_filename;
if ($filter_out) {
   $filtered_dosages_filename = $dosages_filename . "_filtered";
  $duplicatesearch_command .= " -filtered_out $filtered_dosages_filename";
print STDERR "$filtered_dosages_filename \n";
}
#exit;

print "$duplicatesearch_command \n";
###############################################################################


####   Run duplicatesearch :   ################################################
system "$duplicatesearch_command";
###############################################################################

####   Run clusterer.pl    #####################################################

my $cluster_command = $bindir . "/clusterer.pl " . " -in duplicatesearch.out -out cluster.out ";
$cluster_command .= " -cluster $cluster_max_distance  -nofull ";
print $cluster_command, "\n";
#exit();
system $cluster_command;
###############################################################################


####  Store ids of cluster accessions in    ###################################
open my $fh_clusters, "<", "cluster.out";
my %clusterids = ();

while (my $line = <$fh_clusters>) { # each line is one cluster
  next if($line =~ /^\s*#/);
  my @cols = split(" ", $line);
  my $cluster_size = shift @cols;
  my $min_d = shift @cols;
  my $max_d = shift @cols;
  my $n_bad = shift @cols;
  my $rep_id = $cols[0];     # id of the representative of the cluster
  for my $a_cluster_id (@cols) {
    $clusterids{$a_cluster_id} = 1; # all accessions in clusters get stored in @clusterids
  }
}
close $fh_clusters;
###############################################################################


####   Store dosages  #########################################################
print STDERR "before storing dosages\n";

# store individual lines of dosage file in hash
$input_dosages_filename = $filtered_dosages_filename if($filter_out);
open my $fh_dosages, "<", "$input_dosages_filename";
open my $fhout, ">", "$output_dosages_filename" or die "Couldn't open $output_dosages_filename for writing.\n";
open my $fhpheno_out, ">", 'avg_' . $phenotype_filename if(defined $phenotype_filename);
# store ids and genotypes of clusters (size >= 2), and output singletons.
my %id_gts  = ();  # key: ids; value: array ref of dosages (0,1,2,NA) 
my $first_line = <$fh_dosages>;
print $fhout $first_line;
while (my $line = <$fh_dosages>) {
  if ($line =~ /^\s*(#|MARKER)/) {
    print $fhout $line;
    next;
  }
  my @cols = split(" ", $line);
  my $id = shift @cols;
  if (exists $clusterids{$id}) {
    $id_gts{$id} = \@cols;
  } else {
    # print STDERR "ABCD: $id \n";
    print $fhout $line;
    if(defined $phenotype_filename){
    print $fhpheno_out "$id $id  ", $accid_phenotype{$id}, "\n"; # output the phenotype of a non-duplicated accession
  }
  }
}
close($fh_dosages);
print STDERR "after storing dosages\n";


###############################################################################
# my $file_delete_success = unlink($cleaned_dosages_filename);
# warn "Deleting of $cleaned_dosages_filename failed.\n" if($file_delete_success != 1);

my @duplicate_lines = ();
####   Cluster members vote on correct genotypes   ############################
open  $fh_clusters, "<", "cluster.out";
if ($vote  and  $cluster_fraction == 0) { # cluster members vote, and then output just representative id with 'elected' dosages.
  print STDERR "VOTE.\n";
} else {
  print STDERR "DON'T VOTE.\n";
}
while (my $line = <$fh_clusters>) { # each line is one cluster
  next if($line =~ /^\s*#/);

  my @cols = split(" ", $line);
  next if(scalar @cols  < 11);
  my $cluster_size = shift @cols;

  my $min_d = shift @cols;
  my $avg_d = shift @cols;
  my $max_d = shift @cols;
  my $min_intraextra_d = shift @cols;

  my $n_near1 = shift @cols;
  my $n_near2 = shift @cols;
  my $n_big_intra_d = shift @cols;
  my $n_missing_dist = shift @cols;

  my $rep_id = $cols[0];
  die if(scalar @cols != $cluster_size);
  if(exists $accid_phenotype{$rep_id}){
    my $cluster_phenotype = sum(map($accid_phenotype{$_}, @cols))/$cluster_size;
    print $fhpheno_out "$rep_id $rep_id  $cluster_phenotype\n";
  }

  if ($vote  and  $cluster_fraction == 0) {
    my $rep_id = $cols[0];   # id of the representative of the cluster
    my $cluster_id = $rep_id; # . "_size" . $cluster_size; 
    my $elected_gts = vote(\@cols, \%id_gts);
    print $fhout "$cluster_id  ", join(" ", @$elected_gts), "\n";
    
  } else { # don't vote - just use the first one
    my $rep_id = shift @cols; # id of the representative of the cluster
    my $cluster_id = $rep_id; # . "_size" . $cluster_size;
    print $fhout "$cluster_id  ", join(" ", @{$id_gts{$rep_id}}), "\n"; # output representative and its dosages.
    
    for my $an_id (@cols) {
      my $dupe_acc_line = 'DDDD_' . "$an_id  " . join(" ", @{$id_gts{$an_id}}) . "\n"; # other cluster members and dosages.
      push @duplicate_lines, $dupe_acc_line;
    }
  }
}
close $fh_clusters;
###############################################################################
# output some fraction of duplicates:
@duplicate_lines = shuffle @duplicate_lines;
my $n_duplicates_to_output = int($cluster_fraction * scalar @duplicate_lines + 0.5);
for my $i (1..$n_duplicates_to_output) {
  print $fhout $duplicate_lines[$i];
}
close $fhout;


# while (my $line = <$fh_clusters>) { # each line is one cluster
#   next if($line =~ /^\s*#/);
#   my @cols = split(" ", $line);
#   my $min_md_id = $cols[4];
#   my $min_md_count = $id_mdcount{$min_md_id};
#  # if (1) { # uniquified file has accession with the least missing data.
#     for my $a_cluster_id (@cols[4..$#cols]) {
#       $clusterids{$a_cluster_id} = 1; # all accessions in clusters get stored in @clusterids
#       if (0) {
# 	if ($id_mdcount{$a_cluster_id} < $min_md_count) { # find the member of cluster with the least missing data
# 	  $min_md_count =  $id_mdcount{$a_cluster_id};
# 	  $min_md_id = $a_cluster_id;
# 	}
#       }
#     }
#     print $id_line{$min_md_id};
#   # } else { # uniquified file just has first-occurring accession in cluster 
#   #   print $id_line{$cols[4]};
#   # }
# }

# print the lines of the singletons (accessions in clusters of size 1)
# my @sorted_ids = sort {$a cmp $b} keys %id_gts; # sort keys and output in this order for repeatability
# #while (my($an_id, $a_line) = each %id_line) {
# for my $an_id (@sorted_ids) {
#   if (! exists $clusterids{$an_id}) {
#     print "$an_id  ", join(" ", $id_gts{$an_id}), "\n";
#   }
# }

# ************************************************

sub vote{
  my $cluster_ids = shift; # array ref holding the ids of the accessions in the cluster.
  my $id_dosages = shift; # hash ref. key: id, value: array ref of dosages for this id.
  # print STDERR "in vote. new cluster:  \n";
  my $first_id = $cluster_ids->[0];
  my $first_dosages = $id_dosages->{$first_id};
  my $n_markers = scalar @$first_dosages;
  my @elected_dosages = (0) x $n_markers;
  my $cluster_size = scalar @$cluster_ids;

  # my $md_count_novote_2 = 0;
  # my $md_count_vote_2 = 0;
  # my $md_count_novote_gt2 = 0;
  # my $md_count_vote_gt2 = 0;
  
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
	#  $md_count_novote_gt2++; # count the missing gts in the cluster
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
    #  $md_count_vote_gt2++ if($e eq $missing_str);
      $elected_dosages[$i] = $e;
      #   print STDERR "$cluster_size    $i  ", join(" ", @$mv), "    $e \n" if($cluster_size > 20  and  $mv->[3] > 0);
    }
     my $X_count = sum(map( (($_ eq 'X')? 1 : 0), @elected_dosages)); #  =~ tr/X/X/;
   # print STDERR "$cluster_size  $first_id   $md_count_novote_gt2    ", $md_count_novote_gt2/$cluster_size, "  $md_count_vote_gt2   $X_count\n";
  } else {			# cluster size == 2
  #  printf(STDERR "%30s  %s\n", $first_id, join(" ", @{$first_dosages}[0..63]));
    my $second_id = $cluster_ids->[1];
  #  printf(STDERR "%30s  %s\n", $second_id, join(" ", @{$id_dosages->{$second_id}}[0..63]));
 #   my ($okagree_count, $oneX_count, $XX_count, $disagree_count) = (0, 0, 0, 0);
    while (my($i, $d1) = each @$first_dosages) {

      my $d2 = $id_dosages->{$second_id}->[$i];
  #    $md_count_novote_2++ if($d1 eq $missing_str);
  #    $md_count_novote_2++ if($d2 eq $missing_str);
      if ($d1 eq $d2) {
	$elected_dosages[$i] = $d1;
#	if($d1 eq $missing_str){
# 	  $XX_count++;
# #	  $md_count_vote_2++;
# 	}else{
# 	  $okagree_count++;
# 	}
      } elsif ($d1 eq $missing_str) { # d1 is missing, use d2 (which is not missing)
	$elected_dosages[$i] = $d2;
	#  $oneX_count++;
      } elsif ($d2 eq $missing_str) {	# d2 is missing, use d1 (which is not missing)
#	$oneX_count++;
	$elected_dosages[$i] = $d1;
      } else {		    # neither is missing, and they disagree; use $missing_str
	$elected_dosages[$i] = $missing_str;
#	$md_count_vote_2++;
#	$disagree_count++;
      }
    }
  #  printf(STDERR "%30s  %s\n", $first_id . '_grp', join(" ", @elected_dosages[0..63]));
  # my $X_count = sum(map( (($_ eq 'X')? 1 : 0), @elected_dosages)); #  =~ tr/X/X/;
  #   my $nv_Xcount = $oneX_count + 2*$XX_count;
  #   my $v_Xcount = $XX_count + $disagree_count;
  #   print STDERR "$cluster_size  $first_id   $md_count_novote_2    ", $md_count_novote_2/$cluster_size, "  $md_count_vote_2   $X_count   ", scalar @$first_dosages,  "  $okagree_count  $oneX_count  $XX_count  $disagree_count  $nv_Xcount $v_Xcount\n";
  } # end cluster size == 2 branch
  return \@elected_dosages;
}

sub usage_message{
  print "Usage: uniquify -i <input filename> [-o <output filename>] [-m <max allowed fraction marker missing data>].\n";
}

