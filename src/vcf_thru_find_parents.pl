#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename 'dirname';
use List::Util qw 'min max sum';
use File::Spec 'splitpath';

# start with a vcf file -> dsgm file,
# then run duplicate_search
# then clusterer
# then uniquify
# then find_parents
#    either without pedigrees
#    or if there is a pedigree file, run find_parents in mode 1 (-alt 1)

my $vcf_path = undef;
my $pedigree_path = undef;
my $dcluster = undef; # max agmr for which edges will be made in clustering graph.

GetOptions(
	   'vcf_in=s' => \$vcf_path,
	   'pedigree_filename|pedigree_path=s' => \$pedigree_path,
	   'cluster_distance|dcluster=f' => \$dcluster,	   
	  );

my ($v, $vcfdir, $vcf_filename) = File::Spec->splitpath($vcf_path);
print "$v  $vcfdir  $vcf_filename\n";

my ($pv, $peddir, $pedigree_filename) = (defined $pedigree_path)? File::Spec->splitpath($pedigree_path) : (undef, undef, undef);

my $dsgm_filename = $vcf_filename . '.dsgm'; 
my $vcf2dsgm_command = "vcf_to_dsgm  -vcf $vcf_path  -out  $dsgm_filename "; # just use mostly defaults for now.
print "running:  $vcf2dsgm_command \n";
system $vcf2dsgm_command;


my $dsout_filename = $dsgm_filename . ".dsout";
my $ds_command = "duplicate_search  -in $dsgm_filename  -out $dsout_filename  "; # mostly defaults
print "running: $ds_command \n";
system $ds_command;


my $clstrout_filename = $dsout_filename . ".clstrout";
my $clstr_command = "clusterer  -in $dsout_filename  -out $clstrout_filename "; # mostly defaults
$clstr_command .= "  -dclust $dcluster " if(defined $dcluster);
print "running:  $clstr_command \n";
system "$clstr_command";


my $udsgm_filename = 'u_' . $dsgm_filename;
my $upedigree_filename = (defined $pedigree_filename)? 'u_' . $pedigree_filename : undef;
my $uniq_command = "uniquify  -dosages_in $dsgm_filename  -dosages_out $udsgm_filename  -cluster_filename $clstrout_filename ";
$uniq_command .= "  -pedigrees_in $pedigree_path  -pedigrees_out $upedigree_filename " if(defined $pedigree_path);
print "running: $uniq_command \n";
system "$uniq_command";


my $fpout_filename = "$udsgm_filename" . '.fpout';
my $find_parents_command = "find_parents  -in $udsgm_filename  -out $fpout_filename ";  
if(defined $upedigree_filename){
   $find_parents_command .= " -alt 1 -ped $upedigree_filename \n";
 }
print "running:  $find_parents_command \n";
system "$find_parents_command \n";


