package Plot;
use strict;
use warnings;
use Moose;
#use Mouse;
use namespace::autoclean;
use Carp;
use Scalar::Util qw (looks_like_number );
use List::Util qw ( min max sum );
use POSIX qw ( floor ceil );

has histograms => (
		   isa => 'Object',
		   is => 'rw',
		   required => 1,
		  );

has width => (
	      isa => 'Int',
	      is => 'rw',
	      required => 0,
	      default => 640,
	     );

has height => (
	       isa => 'Int',
	       is => 'rw',
	       required => 0,
	       default => 480,
	      );

has key_horiz_position => (
			   isa => 'Str',
			   is => 'rw',
			   required => 0,
			   default => 'center',
			  );

has key_vert_position => (
			   isa => 'Str',
			   is => 'rw',
			   required => 0,
			   default => 'top',
			  );

has ymin => (
	     isa => 'Maybe[Num]',
	     is => 'rw',
	     required => 0,
	     default => undef,
	    );

has ymax => (
	     isa => 'Maybe[Num]',
	     is => 'rw',
	     required => 0,
	     default => undef,
	    );

has ymin_log => (
		 isa => 'Maybe[Num]',
		 is => 'rw',
		 required => 0,
		 default => undef,
		);

has ymax_log => (
		 isa => 'Maybe[Num]',
		 is => 'rw',
		 required => 0,
		 default => undef,
		);

has line_width => (
		   isa => 'Num',
		   is => 'rw',
		   required => 0,
		   default => 1,
		  );

has relative_frame_thickness => ( # thickness of line framing the plot relative to histogram linewidth
				 isa => 'Num',
				 is => 'rw',
				 required => 0,
				 default => 1.5,
				);

has x_axis_label => (
		     isa => 'Maybe[Str]',
		     is => 'rw',
		     required => 0,
		     default => undef,
		    );

has y_axis_label => (
		     isa => 'Maybe[Str]',
		     is => 'rw',
		     required => 0,
		     default => undef,
		    );

has log_y => (
	      isa => 'Bool',
	      is => 'rw',
	      required => 0,
	      default => 0,
	     );

has vline_position => (
		       isa => 'Maybe[Num]',
		       is => 'rw',
		       required => 0,
		       default => undef,
		      );

has output_filename => (
			isa => 'Str',
			is => 'rw',
			required => 0,
			default => 'histogram.png',
		       );

__PACKAGE__->meta->make_immutable;

1;
