package Chromosome;		# phased genotype info for 1 chromosome
use strict;
use warnings;
#use Moose;
use Mouse;
use namespace::autoclean;
use Carp;
use List::Util qw ( min max sum );
#use Marker;

has i_chrom => (
		   isa => 'Int',
		   is => 'ro',
		   required => 1,
	       );

has genotypes => (
		   isa => 'ArrayRef[Maybe[Int]]',
		   is => 'rw',
		   required => 1,
		 );

sub add_genotype{
  my $self = shift;
  my $gt = shift;
  push @{$self->genotypes()}, $gt;
}

sub printxxx{
  my $self = shift;
  my ($n00, $n01, $n10, $n11, $nX, $nother) = (0, 0, 0, 0, 0, 0);
  my $gts = $self->genotypes();
  for my $agt (@$gts){
    if($agt eq 'X'){
      $nX++;
    }elsif($agt == 0){
      $n00++;
    }elsif($agt == 1){
      $n01++;
    }elsif($agt == -1){
      $n10++;
    }elsif($agt == 2){
      $n11++;
    }else{
      $nother++;
    }
  }
  printf(STDOUT "chrom: %5d   n00 etc:  %5d   %5d %5d   %5d   %5d %5d\n", $self->i_chrom(), $n00, $n01, $n10, $n11, $nX, $nother);
}
1;


