#!/usr/bin/perl -w
use strict;
use List::Util qw(min max sum);
use Getopt::Long;


my $input_filename = undef;
my $output_filename = 'bt.out';
my $the_col = 18;	# unit-based
my $verbose = 0;
my $max_solns_out = 2;

GetOptions(
	   'input_file=s' => \$input_filename, 
	   'output_file=s' => \$output_filename,
	   'column=i' => \$the_col, # 
	   'verbose!' => \$verbose,
	   'solutions=i' => \$max_solns_out,
	  );
$the_col--; # now zero-based

open my $fhin, "<", "$input_filename" or die "Couldn't open $input_filename for reading.\n";
my %accid_solns = ();
while (my $line = <$fhin>) {
  next if($line =~ /^\s*#/);
   my @cols = split(" ", $line);
   my $value = $cols[$the_col];
  my $accid = shift @cols;
  my $parents = $cols[0] . " " . $cols[1];
  my $cross_type = ($cols[0] eq $cols[1])? "self" : "bip";
  my $x = join(" ", @cols);
  my $vc = $value . "\t" .  "$cross_type $x ";
  if (exists $accid_solns{$accid}) {
    push @{$accid_solns{$accid}}, $vc; # "$parents  $value"; 
  } else {
    $accid_solns{$accid} = [$vc];
  }
}
close $fhin;

open my $fhout, ">", "$output_filename" or die "Couldn't open $output_filename for writing.\n";
for my $anid (keys %accid_solns) {
  my $solutions = $accid_solns{$anid};
  my @sorted_solns = sort {val($a) <=> val($b)} @$solutions;
  my $n_solns_out = min(scalar @sorted_solns, $max_solns_out);
  print $fhout "$anid  ", join("  ", @sorted_solns[0..$n_solns_out-1]), "\n";
 # print "  ", $sorted_solns[1] if(scalar @sorted_solns > 1);
 # print "\n"; # join("  ", @sorted_solns), "\n";
}


sub val{
  my $vc = shift;
  my ($value, $xx) = split("\t", $vc);
  return $value;
}
