package Hdata;
use strict;
use warnings;
use Moose;
#use Mouse;
use namespace::autoclean;
use Carp;
use Scalar::Util qw (looks_like_number );
use List::Util qw ( min max sum );
use POSIX qw ( floor ceil );

# to hold the data for one histogram (typically from one column, but could depend on >1 col),
# together with mean, stddev, etc.

# data values, sorted but not binnedx

has data_array => (
                  isa => 'ArrayRef',
                  is => 'rw',
                  default => sub { [] },
		  );

has log0_count => ( # if histogramming log of col values, this counts occurences of value=0
		   # these will be added to the underflow bin
		   isa => 'Num',
		   is => 'rw',
		   default => 0
		   );

######  summary statistics (mean, stddev, etc.) #########################

has n_points => (
                 isa => 'Int',
                 is => 'rw',
                 default => 0,
    );

has n_undefined => (
    isa => 'Int',
    is => 'rw',
    default => 0,
    );

has min => (
             isa => 'Maybe[Num]',
             is => 'rw',
             default => undef,
            );

has max => (
             isa => 'Maybe[Num]',
             is => 'rw',
             default => undef,
            );

# has range => (                  # min and max numbers in
#               isa => 'ArrayRef[Maybe[Num]]',
#               is => 'rw',
#               default => sub { [undef, undef] },
#              );

has sum => (
             isa => 'Maybe[Num]',
             is => 'rw',
             default => 0,
            );

has sumsqr => (
             isa => 'Maybe[Num]',
             is => 'rw',
             default => 0    ,
            );

has mean => (
             isa => 'Maybe[Num]',
             is => 'rw',
             default => undef,
            );
has stddev => (
               isa => 'Maybe[Num]',
               is => 'rw',
               default => undef,
              );
has stderr => (
               isa => 'Maybe[Num]',
               is => 'rw',
               default => undef,
              );
has median => (
               isa => 'Maybe[Num]',
               is => 'rw',
               default => undef,
    );

has label => (
    isa => 'Maybe[Str]',
    is => 'rw',
    default => undef,
    );

#########################################################################

######  binned data #####################################################

has bin_counts => (
                   isa => 'ArrayRef',
                   is => 'rw',
                  );

has bin_centers => (
                    isa => 'ArrayRef',
                    is => 'rw',
                   );

has underflow_count => (
                        isa => 'Num',
                        is => 'rw',
                       );
has overflow_count => (
                       isa => 'Num',
                       is => 'rw',
                      );

sub add_value{
   my $self = shift;
   my $value = shift;
   if(defined $value){
   push @{$self->data_array()}, $value;
   $self->{sum} += $value;
   $self->{sumsqr} += $value*$value;
   }else{
       $self->{n_undefined}++;
   }
}

sub sort_etc{
   my $self = shift;

   my $n_points =  scalar @{ $self->data_array() };
   if ($n_points > 0) {
      $self->{data_array} = [sort {$a <=> $b} @{$self->{data_array}}];
      $self->min( $self->data_array()->[0] );
      $self->max( $self->data_array()->[-1] );
  
      $self->mean( $self->sum()/$n_points );
      my $mean = $self->sum()/$n_points;
      my $variance = $self->sumsqr()/$n_points - $mean*$mean;
      my $stddev = sqrt($variance);
      my $stderr = $stddev/sqrt($n_points);
      $self->n_points( $n_points );
      $self->mean( $mean );
      $self->stddev( $stddev );
      $self->stderr( $stderr );

      $self->median( $self->calculate_median() );
      # if ($n_points % 2 == 0) {
      #    my $mid = int($n_points/2);
      #    $self->median( 0.5*($self->{data_array}->[$mid] + $self->{data_array}->[$mid+1]) );
      # } else {
      #    $self->median( $self->{data_array}->[ int($n_points/2) ] );
      # }

   }
 }

sub add_to_log0_count{
  my $self = shift;
  my $increment = shift // 1;
  $self->{log0_count} += $increment;
}

sub calculate_median{
  my $self = shift;
  # if histogramming log of values, include log(0) pts in count for median:
  my $median_n_points = $self->n_points() + $self->log0_count();
  my $median_index;
  my $median = undef;
  if($median_n_points % 2 == 0){ # even
    $median_index = $median_n_points/2 - $self->log0_count();
    $median = 0.5*($self->{data_array}->[$median_index] + $self->{data_array}->[$median_index-1]) if($median_index >= 0);
  }else{
    my $median_index = ($median_n_points-1)/2 - $self->log0_count();
    $median = $self->{data_array}->[$median_index] if($median_index >= 0);
  }
  return $median;
}

1;
