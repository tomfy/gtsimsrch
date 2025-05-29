#!/usr/bin/perl -w
use strict;
use Getopt::Long;

# read in a cluster file (output of agmr_cluster ) and a
# pedigree file (whitespace separated, accession, Fparent, Mparent ids in last 3 cols of each line)
# output to stdout a pedigree file with the parental ids replaced with the representative of the cluster

my $cluster_file = undef;
my $pedigree_file = undef;
my $output_file = undef;
my $progeny_col = -3;
my $Fparent_col = -2;
my $Mparent_col = -1;
my $nskip = 0; # skip this many lines at top of pedigree file.

my $repid_column = 10; # unit-based, the column with the id of the representative (min missing data) accession of duplicate group.
$repid_column--;       # now zero-based

GetOptions(
	   'cluster_file=s' => \$cluster_file,
	   'pedigree_file=s' => \$pedigree_file,
	   'output_file=s' => \$output_file,
	   'progeny_column=i' => \$progeny_col,
	   'Fparent_column=i' => \$Fparent_col,
	   'Mparent_column=i' => \$Mparent_col,
	   'skip|nskip=i' => \$nskip,
	  );

if(!defined $cluster_file  or  !defined $pedigree_file){
if(!defined $cluster_file){
print STDERR "Cluster file is undefined; it must be specified\n";
}
if(!defined $pedigree_file){
print STDERR "Pedigree file is undefined; it must be specified\n";
}
exit;
}
if(!defined $output_file){
  $output_file = 'u_' . $pedigree_file;
}

# want col 1 to be leftmost, -1 rightmost
$progeny_col-- if($progeny_col > 0);
$Fparent_col-- if($Fparent_col > 0);
$Mparent_col-- if($Mparent_col > 0);

###  store cluster info in hash  ###
my %clusterid_repid = (); # key: an id in a cluster, value: the representative id for that cluster
my %repids = (); # keys are ids of the representative ids of the clusters (with size >= 2)
open my $fhin, "<", "$cluster_file";
while (my $line = <$fhin>) {
  next if($line =~ /^\s*#/);
  next if($line =~ /^\s*$/);
  # print STDERR "[ $line ]";
  my @cols = split(" ", $line);
  # print STDERR $cols[0], "\n";
  my $the_repid = $cols[$repid_column];
  $repids{$the_repid} = 1;
  for my $id (@cols[$repid_column..$#cols]) {
    # print STDERR "$id $the_repid\n";
    $clusterid_repid{$id} = $the_repid; # key is id of a cluster member, and value is the id of the representative member of group of duplicates.
  }
}
close $fhin;
print STDERR "# n duplicate groups: ", scalar keys %repids, "  n ids in duplicate groups: ", scalar keys %clusterid_repid, "\n";

###  now read pedigree file and replace each accession id with the representative of the duplicate group to which it belongs  ###
open $fhin, "<", "$pedigree_file";
open my $fhout, ">", "$output_file";
for(1..$nskip){ <$fhin>; } # skip first nskip lines.

my $pedigrees_in_count = 0;
my $progeny_NA_count = 0;
my $duplicate_progeny_count = 0;
my $pedigrees_out_count = 0;

my %uniqprogeny = ();
while (my $line = <$fhin>) {
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
      print $fhout "$progeny_id  $Fpar_id  $Mpar_id\n";
    }
  } else {
    $progeny_NA_count++;
  }
}
close $fhin;
print STDERR "# $duplicate_progeny_count \n";

print STDERR "# pedigrees read in: $pedigrees_in_count\n";
print STDERR "# progeny id is NA count: $progeny_NA_count\n";
print STDERR "# duplicate progeny count: $duplicate_progeny_count\n";
print STDERR "# pedigrees output: $pedigrees_out_count\n";
