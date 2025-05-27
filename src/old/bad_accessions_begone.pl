#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw(min max sum);

# read a dosage matrix file,
# remove accessions which have excessive (>50% by default) missing data,
# and output with same format as input file.

my $input_dosages_filename = undef;
my $output_dosages_filename = undef;
my $max_acc_missing_data_fraction = 0.5;

GetOptions(
	   'input_file=s' => \$input_dosages_filename,
	   'output_file=s' => \$output_dosages_filename,
	   'max_md_fraction=f' => \$max_acc_missing_data_fraction,
	  );

if (!defined $input_dosages_filename) {
  print "Must supply input file name.\n";
  usage_message();
  exit;
}

print STDERR "bad_accessions_begone output filename $output_dosages_filename \n";

open my $fhin, "<", "$input_dosages_filename" or die "Couldn't open $input_dosages_filename for reading.\n";
$output_dosages_filename = $input_dosages_filename . "_bad_accessions_removed" if(!defined $output_dosages_filename);


####   remove accessions with excessive missing data   #########################
open my $fhout, ">", "$output_dosages_filename";
my $n_markers = undef;
my $n_bad_accessions = 0;

printf(STDOUT "# Removing accessions with > %5.2f percent missing data.\n",
       $max_acc_missing_data_fraction*100);
while (my $line_in = <$fhin>) {
  if ($line_in =~ /^\s*#/) {
    print $fhout $line_in; # comment lines just pass through unaffected
  } elsif ($line_in =~ /^MARKER/) {
    print $fhout $line_in;
    my @markers = split(" ", $line_in);
    $n_markers = scalar @markers  - 1;
  } else {
    die "File lacks line with MARKER and marker ids\n" if(! defined $n_markers);
    my $n_bad = () = $line_in =~ /\sNA/g; # count markers with missing data
    if ($n_bad <= $max_acc_missing_data_fraction*$n_markers) {
      print $fhout $line_in;	# output a good accession
    } else {
      $n_bad_accessions++;
    }
  }
}
print STDOUT "# $n_bad_accessions accessions removed due to excessive missing data.\n";
#printf $fhout "# $n_bad_accessions accessions removed due to > %5.2f percent missing data\n",
#  $max_acc_missing_data_fraction*100;
close $fhout;

sub usage_message{
  print "Usage: bad_accessions_begone.pl -i <input filename> [-o <output filename>] [-m <max allowed fraction missing data>].\n";
}
