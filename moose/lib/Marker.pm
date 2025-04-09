package Marker;		# phased genotype info for 1 diploid marker
use strict;
use warnings;
#use Moose;
use Mouse;
use namespace::autoclean;
use Carp;
use List::Util qw ( min max sum );

# has i_chrom => (
# 		   isa => 'Int', 
# 		   is => 'ro',
# 		   required => 1,
# 		  );
# has id => (
# 		   isa => 'Str',
# 		   is => 'ro',
# 		   required => 1,
# 		  );
# has position => (
# 		   isa => 'Int',
# 		   is => 'ro',
# 		   required => 1,
# 		  );

has pgt => (
		     isa => 'Int', # 00->0, 01->1, 10->2, 11->3, missing -> 4
		     is => 'ro',
		     required => 1,
		    );

1;
