#!/usr/bin/perl -w
use strict;

my $partial_filename = shift; # plnk -> plnk.mdist plnk.mdist.id

my $dist_filename = $partial_filename . ".mdist";
my $id_filename = $partial_filename . ".mdist.id";

my $filename_out = shift;
my $max_distance = shift // 1.1;

# my $filename_out = $partial_filename . ".dists";

#read in ids
my @ids = ();
my $badline_count = 0;
open my $fhid, "<", "$id_filename";
while(<$fhid>){
  next if(/^\s*#/);
  my $id;
if(/^\s*(\S+)/){
  $id = $1;
  push @ids, $id;
}else{
  $badline_count++;
  print "empty line in input file $id_filename.\n";
}

}
print STDERR "# read ", scalar @ids, " accession ids.\n";
close $fhid;

open my $fhdist, "<", "$dist_filename";
open my $fhout, ">", "$filename_out";
my $line_number = 0;
while(my $line = <$fhdist>){
  # next if($line =~ /^\s*#/);
  my $id1 = $ids[$line_number];
  my @distances = split(" ", $line);
  for(my $i = 0; $i < $line_number; $i++){
      my $distance = $distances[$i];
      if($distance <= $max_distance){
#	my $id2 = $ids[$i];
	printf($fhout "%s  %s  %10.8f \n", $id1, $ids[$i], $distance);
      }
    }
  $line_number++;
}
close $fhdist;
close $fhout;
