package Chromosome;		# phased genotype info for 1 chromosome
use strict;
use warnings;
#use Moose;
use Mouse;
use namespace::autoclean;
use Carp;
use List::Util qw ( min max sum );
use Marker;

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

1;


