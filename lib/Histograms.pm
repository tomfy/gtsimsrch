package Histograms;
use strict;
use warnings;
use Moose;
#use Mouse;
use namespace::autoclean;
use Carp;
use Scalar::Util qw (looks_like_number );
use List::Util qw ( min max sum );
use POSIX qw ( floor ceil );
use Hdata;

my $x_notight_factor = 2.0;
my $x_tight_factor = 0.05;

use constant  BINWIDTHS => {
			    100 => [80, 125],
			    125 => [100, 200],
			    200 => [125, 250],
			    250 => [200, 400],
			    400 => [250, 500],
			    500 => [400, 800],
			    800 => [500, 1000]
			   };

#######  which file and column of file to read  #########################


has data_fcol => ( # string representing files and column(s) (unit based) in which to find the numbers to be histogrammed.
                  isa => 'Str', # e.g. '0v1:2,3,4; 0v2:4,7,10' -> histogram cols 2,3,4 of file 0v1, and cols 4,7,10 of file 0v2
                  is => 'rw',
                  required => 1,
                 );

has filecol_specifiers => ( # one entry specifying file and column (e.g. 'x.out:3' for each histogram.
			   # also e.g. file1:3+4  or file1:3-5"difference of cols 3,5"
                           isa => 'ArrayRef',
                           is => 'rw',
                           required => 0,
			  );

has histograms_to_plot => ( # if you request n histograms, this will default to [1,1,1,...1] (i.e. n 1's]
			   # but if one of these is set to 0, then corresponding histogram will not be plotted
			   isa => 'ArrayRef',
			   is => 'rw',
			   required => 0,
			  );

# has data_files_and_columns => ( # Hashref. keys are files, values strings specifying cols to be histogrammed.
#                                isa => 'HashRef',
#                                is => 'rw',
#                                required => 0,
#                               );

#########################################################################


#######  the raw (unbinned) data  #######################################

has data_type => (
                  isa => 'Maybe[Str]', # 'integer' or 'float'
                  is => 'rw',
                  default => undef,
                 );

has filecol_hdata => (
                      # a hashref of Hdata objects, one for each histogram,
                      # keys are strings typically with filenames and column numbers, e.g. 'x.out:1' 
                      # but potentially something like 'x.out:3+5' (histogram of the sum of col 3 and col 5 values)
                      # plus one for all of them pooled, with key 'pooled'
                      isa => 'HashRef',
                      is => 'rw',
                      default => sub { {} },
                     );

# has data_hash => ( # needed?
#              isa => 'HashRef',
#              is => 'ro',
#              default => sub { {} },
#             );
#########################################################################

######  binning specifiers ##############################################
has tight => (
		   isa => 'Bool',
		   is => 'rw',
		   required => 0,
		   default => 1,
		   );

has binwidth => (
                 isa => 'Maybe[Num]',
                 is => 'rw',
                 required => 0,
		 default => undef,
                );

has lo_limit => (               # low edge of lowest bin
                 isa => 'Maybe[Num]',
                 is => 'rw',
                 required => 0,
		 default => undef,
                );

has hi_limit => (               # high edge of highest bin
                 isa => 'Maybe[Num]',
                 is => 'rw',
                 required => 0,
		 default => undef,
                );

has n_bins => (
               isa => 'Maybe[Num]',
               is => 'rw',
               required => 0,
	       default => undef,
              );

has titles => (
	       isa => 'Maybe[ArrayRef[Str]]',
	       is => 'rw',
	       required => 0,
	       default => undef,
	      );

has max_bin_y => (
		  isa => 'Num',
		  is => 'rw',
		  required => 0,
		  default => 0,
		  );

#########################################################################


#########################################################################

# around BUILDARGS => sub{
#
# };
# don't forget the ';' here!

sub BUILD{
  my $self = shift;

  $self->set_filecol_specs(); # construct filecol_specifiers (e.g. ['x.out:3', 'x.out:5'] from 'x.out:3,5'
  my @histos_to_plot = ((1) x scalar @{$self->filecol_specifiers()});
  $self->histograms_to_plot(\@histos_to_plot); # default to plotting all (but not pooled)
  $self->load_data_from_file();
  # print "lo_limit, hi_limit, binwidth: ",  $self->lo_limit() // 'undef', " ", $self->hi_limit() // 'undef', " ", $self->binwidth() // 'undef', "\n";

  #  if (!(defined $self->lo_limit()  and  defined $self->hi_limit()  and  defined  $self->binwidth())) 
  # if (1 or !defined $self->binwidth()) {
  if(!(defined $self->lo_limit()  and  defined $self->hi_limit()  and  defined  $self->binwidth())){
    my ($datamin, $datamax) = $self->data_min_max();
    print "$datamin $datamax  ", $self->tight(), "\n";
    if($self->tight()){
      if(!defined $self->lo_limit()){
      if($datamin >=0){
	$self->lo_limit(0);
      }else{
	$self->lo_limit($datamin - $x_tight_factor*($datamax-$datamin));
      }
    }
      $self->hi_limit($datamax + $x_tight_factor*($datamax-$datamin)) if(!defined $self->hi_limit());;
    }else{
      if(!defined $self->lo_limit()  and  $datamin >= 0){
	$self->lo_limit(0);
      }
    }
    $self->auto_bin($self->lo_limit(), $self->hi_limit(), $self->binwidth());
  }

  my $n_bins = int( ($self->hi_limit() - $self->lo_limit())/$self->binwidth() ) + 1;
  $self->n_bins($n_bins);

  # print "In BUILD: ", $self->lo_limit(), "  ", $self->hi_limit(), "  ", $self->binwidth(), "  ", $self->n_bins(), "\n";
}

sub load_data_from_file{
  my $self = shift;

  my %filecol_hdata = ('pooled' => Hdata->new());

  my $integer_data = 1;
  for my $histogram_id (@{$self->filecol_specifiers()}) { # $histogram_id e.g. 'file_x:4' or 'file_x:5+6"sum_of_5and6"'
    $filecol_hdata{$histogram_id} = Hdata->new();
    #   print "histogram id: [$histogram_id] \n";
    if ($histogram_id =~ /^([^:]+) [:] (.*) /x) { #
      my ($datafile, $col_specifier) = ($1, $2); # $col_specifier is columns specifier, e.g. 8 or 8+10 or 8u3 = 9-2"x-z"
      my $label = ($col_specifier =~ s/\" (.*) \"//x)? $1 : $histogram_id; # capture label and delete from $col_specifier
      $filecol_hdata{$histogram_id}->label($label);
      $col_specifier =~ s/\s+//g;
      if ($col_specifier =~ /^\d+([+-]\d+)*$/) {
	$integer_data = read_and_sum_signed($datafile, $histogram_id, $col_specifier, \%filecol_hdata);
      } elsif ($col_specifier =~ /^\d+(u\d+)*$/) {
	$integer_data = read_and_union($datafile, $histogram_id, $col_specifier, \%filecol_hdata);
      } elsif ($col_specifier =~ /^\d+[\/]\d+$/) {
	$integer_data = read_and_divide($datafile, $histogram_id, $col_specifier, \%filecol_hdata);
      } elsif ($col_specifier =~ /^\d+[\*]\d+$/) {
	$integer_data = read_and_multiply($datafile, $histogram_id, $col_specifier, \%filecol_hdata);
	  } elsif ($col_specifier =~ /^log\d+$/) {
	$integer_data = read_and_log($datafile, $histogram_id, $col_specifier, \%filecol_hdata);
      } else {
	print STDERR "Unrecognized column and operation specification: $col_specifier. Skipping histogram $histogram_id\n";
      }
    }
  }

  $self->data_type( ($integer_data)? 'integer' : 'float' );
  for my $hdata (values %filecol_hdata) {
    $hdata->sort_etc();
  }
 print "pooled data min, max: ", $filecol_hdata{'pooled'}->min(), "  ", $filecol_hdata{'pooled'}->max(), "\n";
  $self->filecol_hdata( \%filecol_hdata );
}

######## defining the binning #########

sub auto_bin{			# automatically choose binwidth, etc.
  my $self = shift;
  my $lolimit = shift // undef;
  my $hilimit = shift // undef;
  my $binwidth = shift // undef;

  my $pooled_hdata = $self->filecol_hdata()->{'pooled'};
  my @bws = sort( keys %{ BINWIDTHS() } );
  my ($x_lo, $x_hi) = ($pooled_hdata->min(), $pooled_hdata->max()); # the min and max values in the (pooled) data.
  my $n_points = $pooled_hdata->n_points();

  my $i5 = int(0.05*$n_points);
  my $i95 = -1*($i5+1);
  my ($v5, $v95) = ($pooled_hdata->data_array()->[$i5], $pooled_hdata->data_array()->[$i95]);

  my $X = $x_notight_factor;
  my $mid = 0.5*($v5 + $v95);
  my $half_range95 = 0.5*($v95-$v5);

  my ($lo_limit, $hi_limit) = ($mid - $X*$half_range95, $mid + $X*$half_range95);
  $lo_limit = 0 if($x_lo >= 0  and  $lo_limit < 0); # if all data >=0, x scale doesn't include negatives.
  $self->lo_limit( $lolimit // $lo_limit ); # arguments to auto_bin $lolimit and $hilimit override 
  $self->hi_limit( $hilimit // $hi_limit ); # auto calculated values

  if (!defined $binwidth) {
    $binwidth =  $half_range95/$n_points**0.3333;
    ($binwidth, my $bwf) = xyz($binwidth);
    # round binwidth down to nearest allowed value.
    for (my $i = @bws-1; $i >= 0; $i--) {
      my $bw = $bws[$i];
      if ($bw <= $binwidth) {
	$binwidth = $bw;
	last;
      }
    }
    $binwidth *= $bwf;
    if ($self->data_type() eq 'integer') {
      $binwidth = 1 if ($binwidth < 1);
      $binwidth = int($binwidth + 0.5);
    }
  }
  $self->set_binwidth($binwidth);
  print "lo_limit, hi_limit, binwidth from auto_bin: ", $self->lo_limit(), "  ", $self->hi_limit(), "  ", $binwidth,  "\n";
}

sub change_range{
  my $self = shift;
  my $new_lo = shift // undef;
  my $new_hi = shift // undef;
  #   print STDERR "new_lo, new_hi:  ", $new_lo // 'undef', "  ", $new_hi // 'undef', "\n";
  $self->lo_limit($new_lo) if(defined $new_lo);
  $self->hi_limit($new_hi) if(defined $new_hi);
  $self->set_binwidth($self->binwidth());
}

sub expand_range{
  my $self = shift;
  my $factor = shift // 1.2;

  my $mid_x = 0.5*($self->lo_limit() + $self->hi_limit());
  my $hrange = $self->hi_limit() - $mid_x;
  my $lo_limit = $mid_x - $factor*$hrange;
  my $hi_limit = $mid_x + $factor*$hrange;
  $lo_limit = max($lo_limit, 0) if($self->filecol_hdata()->{'pooled'}->min() >= 0); # $self->pooled_hdata()->min
  my $binwidth = $self->binwidth();
  $lo_limit = $binwidth * floor( $lo_limit / $binwidth );
  $hi_limit = $binwidth * ceil( $hi_limit / $binwidth );
  $self->lo_limit( $lo_limit );
  $self->hi_limit( $hi_limit );
  $self->n_bins( int( ($self->hi_limit() - $self->lo_limit())/$self->binwidth() ) + 1 );
  print "after:  ", $self->lo_limit(), '  ', $self->hi_limit(), '  ', $self->binwidth(), '  ', $self->n_bins(), "\n";
}

sub change_binwidth{
  my $self = shift;
  my $n_notches = shift // 1; # 1 -> go to next coarser binning, -1 -> go to next finer binning, etc.
  return if($n_notches == 0);
  my $bw = $self->binwidth();
  if ($n_notches > 0) {
    for (1..$n_notches) {
      my ($bw100, $pow10) = xyz($bw);
      $bw = (BINWIDTHS->{ int($bw100+0.5) }->[1])*$pow10;
    }
    $self->set_binwidth($bw);
  } else {			#  ($n_notches == -1) {
    $n_notches *= -1;
    for (1..$n_notches) {
      my ($bw100, $pow10) = xyz($bw);
      $bw = (BINWIDTHS->{ int($bw100+0.5) }->[0])*$pow10;
    }
    $self->set_binwidth($bw);
  }
}

sub set_binwidth{ # set the binwidth and adjust the lo and hi limits to be multiples of binwidth; set n_bins accordingly.
  my $self = shift;
  my $new_bw = shift;
  my ($lo_limit, $hi_limit) = ($self->lo_limit(), $self->hi_limit());
  $lo_limit = $new_bw * floor( $lo_limit / $new_bw );
  $hi_limit = $new_bw * ceil( $hi_limit / $new_bw );
  if ($self->data_type() eq 'integer') {
    $lo_limit -= 0.5;
    $hi_limit += 0.5;
  }
  my $n_bins =  int( ($self->hi_limit() - $self->lo_limit())/$new_bw ) + 1;
  #  print STDERR "in set binwidth: $lo_limit  $hi_limit  $new_bw $n_bins\n";

  $self->lo_limit($lo_limit);
  $self->hi_limit($hi_limit);
  $self->binwidth($new_bw);
  $self->n_bins($n_bins);
}

sub bin_data{ # populate the bins using existing bin specification (binwidth, etc.)
  my $self = shift;

  my $max_bin_y = -1;
  while (my ($col, $hdata) = each %{$self->filecol_hdata}) {

    my @bin_counts = (0) x ($self->n_bins() + 1);
    my @bin_centers = map( $self->lo_limit() + ($_ + 0.5)*$self->binwidth(), (0 .. $self->n_bins() ) );
    my ($underflow_count, $overflow_count) = (0, 0);
    my ($lo_limit, $hi_limit) = ($self->lo_limit(), $self->hi_limit());

    for my $d (@{$hdata->data_array()}) {
      if ($d < $lo_limit) {
	$underflow_count++;
      } elsif ($d >= $hi_limit) {
	$overflow_count++;
      } else {
	my $bin_number = int( ($d - $lo_limit)/$self->binwidth() );
	#print "$lo_limit  ", $self->binwidth(), "  $d  $bin_number \n";
	$bin_counts[$bin_number]++;
	# $bin_centers[$bin_number] = ($bin_number+0.5)*$self->binwidth()
      }
    }
    $max_bin_y = max(max(@bin_counts), $max_bin_y) if($col ne 'pooled');
    # print "# Col: $col  max_bin_y: $max_bin_y\n";
    my $log0count = $self->filecol_hdata()->{$col}->log0_count();
    $self->filecol_hdata()->{$col}->bin_counts( \@bin_counts );
    $self->filecol_hdata()->{$col}->bin_centers( \@bin_centers );
    $self->filecol_hdata()->{$col}->underflow_count( $underflow_count + $log0count );
    $self->filecol_hdata()->{$col}->overflow_count( $overflow_count );
  }
  $self->max_bin_y($max_bin_y);
}

sub binned_data{
  my $self = shift;
  return ($self->bin_centers(), $self->underflow_counts(), $self->bin_counts(), $self->overflow_counts());
}

sub as_string{
  my $self = shift;
  my $h_string = '';		# the histogram as a string

  #   my @filecol_specs = @{$self->get_filecol_specs()};
  my @filecol_specs = @{$self->filecol_specifiers};
  push @filecol_specs, 'pooled';
  my %fcs_cumulative = ();
  # print STDERR "filecolspecs: ", join("; ", @filecol_specs), "\n";
  my $horiz_line_string .= sprintf("#------------------------------------");
  for my $fcs (@filecol_specs) {
    $horiz_line_string .= "--------------------";
    $fcs_cumulative{$fcs} = 0;
  }
  $horiz_line_string .= "\n";

  $h_string .= sprintf("# data from file:column                " . "%-18s  " x (@filecol_specs) . "\n", @filecol_specs);
  $h_string .= $horiz_line_string;

  $h_string .= sprintf("      n undefined                   ");
  for my $fcspec (@filecol_specs) {
    $h_string .= sprintf("%9.4g           ", $self->filecol_hdata()->{$fcspec}->n_undefined());
  }
  $h_string .= "\n";
  $h_string  .= $horiz_line_string;

  $h_string .= sprintf("     < %6.4g  (underflow)          ", $self->lo_limit());
  for my $fcspec (@filecol_specs) {
    my $val =  $self->filecol_hdata()->{$fcspec}->underflow_count() // 0;
    $fcs_cumulative{$fcspec} += $val;
    $h_string .= sprintf("%9.4g %9.4g ", $val, $fcs_cumulative{$fcspec}); # $self->filecol_hdata()->{$fcspec}->underflow_count() // 0);
  }
  $h_string .= "\n";
  # my $val =  $self->filecol_hdata()->{}->underflow_count() // 0;
  #   $fcs_cumulative{$fcspec} += $val;
  #   $h_string .= sprintf("%9.4g %9.4g", $val, $fcs_cumulative{$fcspec}); # $self->filecol_hdata()->{$fcspec}->underflow_count() // 0);
  # $h_string .= sprintf("%9.4g \n", $self->filecol_hdata()->{pooled}->underflow_count() );

  $h_string .= $horiz_line_string;
  $h_string .= sprintf("# bin     min    center       max     count \n");

  while (my ($i, $bin_center_x) = each @{$self->filecol_hdata()->{pooled}->bin_centers()} ) {
    my $bin_lo_limit = $bin_center_x - 0.5*$self->binwidth();
    my $bin_hi_limit = $bin_center_x + 0.5*$self->binwidth();
    $h_string .= sprintf("    %9.4g %9.4g %9.4g   ",
			 $bin_lo_limit, $bin_center_x,
			 $bin_hi_limit);

    for my $fcspec (@filecol_specs) {
      my $val =  $self->filecol_hdata()->{$fcspec}->bin_counts()->[$i] // 0;
      $fcs_cumulative{$fcspec} += $val;
      $h_string .= sprintf("%9.4g %9.4g ", $val, $fcs_cumulative{$fcspec}); # $self->filecol_hdata()->{$fcspec}->bin_counts()->[$i] // 0);
    }

    $h_string .= "\n"; # sprintf("%9d\n", ($self->filecol_hdata()->{pooled}->bin_counts()->[$i] // 0));
  }
  $h_string .= $horiz_line_string;
  $h_string .= sprintf("     > %6.4g   (overflow)          ", $self->hi_limit());

  for my $fcspec (@filecol_specs) {
    my $val =  $self->filecol_hdata()->{$fcspec}->overflow_count() // 0;
    $fcs_cumulative{$fcspec} += $val;
    $h_string .= sprintf("%9.4g %9.4g ", $val, $fcs_cumulative{$fcspec}); # $self->filecol_hdata()->{$fcspec}->overflow_count() // 0);
  }
  $h_string .= sprintf("\n"); #, $self->filecol_hdata()->{pooled}->overflow_count() );
  $h_string .= $horiz_line_string;

  #  $h_string .= sprintf("# range: [%9.4g,%9.4g]   median: %9.4g\n", @{$self->range()},  $self->median());
  #  for (@col_specs) {
  for my $fcspec (@filecol_specs) {
    my ($f, $cspec) = split(':', $fcspec);
  #  print STDERR "f cspec:  ", $f // 'undef', "  ", $cspec // 'undef', "\n";
    $cspec = '' if(!defined $cspec);
    my @colspecs = split(',', $cspec);
    for my $csp (@colspecs) {
      my $fc = $f . ':' . $csp;
 #     print STDERR "$fcspec  $fc \n";
      my $the_fchd = $self->filecol_hdata()->{$fc};
      $h_string .= sprintf("# file:col %10s  n points: %5d   ", $fc,
			   # $self->filecol_hdata()->{$fc}->n_points());
			   $the_fchd->n_points());
      $h_string .= sprintf("mean: %9.4g   stddev: %9.4g   stderr: %9.4g   median: %9.4g\n",
			   $the_fchd->mean(),
			   $the_fchd->stddev(),
			   $the_fchd->stderr(),
			   $the_fchd->median());
      
    }
  }
  $h_string .= sprintf("# pooled               n points: %5d   ", $self->filecol_hdata()->{pooled}->n_points());
  $h_string .= sprintf("mean: %9.4g   stddev: %9.4g   stderr: %9.4g   median: %9.4g\n",
		       $self->filecol_hdata()->{pooled}->mean(),
		       $self->filecol_hdata()->{pooled}->stddev(),
		       $self->filecol_hdata()->{pooled}->stderr(),
		       $self->filecol_hdata()->{pooled}->median() );
  return $h_string;
}

sub set_filecol_specs{
  my $self = shift;
  my @filecol_specs = ();

  #  my @filecol_specifiers = split(/;/, $self->data_fcol() ); # e.g. '0v1:3,4,5; 0v2:1,5,9' -> ('0v1:3,4,5', '0v2:1,5,9')
  for my $fcs (split(/;/, $self->data_fcol() )) { # split on semi-colon for the different input files
    #   $fcs =~ s/\s+//g; # remove whitespace
    # print STDERR "fcs: ", $fcs, "\n";
    my ($file, $cols) = split(':', $fcs); # split file:spec
    #   my @colspecs = split(',', $cols);
    my @colspecs = split(/[,]+/, $cols); # split on , to get column specs (e.g. 9"hgmr"
    for (@colspecs) {
      #    my $filecolspec = $file . ':' . $_;
      #   print STDERR "filecolspec: [$filecolspec] \n";
      push @filecol_specs, $file . ':' . $_;
    }
  }
  $self->filecol_specifiers( \@filecol_specs ); # example element 'a_filename:8"a_label"', or 'a_filename:10-8'
  #  return \@filecol_specs;
}

#### ordinary subroutines ########

sub xyz{ # express the input number as prod. of 2 factors, one in range [100,1000)
  # the other an int power of 10.
  my $x = shift;
  my $f = 1;
  while ($x >= 1000) {
    $x /= 10;
    $f *= 10;
  }
  while ($x < 100) {
    $x *= 10;
    $f /= 10;
  }
  return ($x, $f);
}



sub read_and_sum_signed{ # signed sum of cols, e.g. 5+8-12  -> add values in cols 5 and 8 and subtract col 12 value.
  my $data_file = shift;
  my $histogram_id = shift;
  my $col_sum_spec = shift;
  my $filecol_hdata = shift;
  open my $fh_in, "<", $data_file or die "Couldn't open $data_file for reading.\n";

  # store columns and signs in arrays
  my @cols_to_sum = ();
  my @signs = ();
  if ($col_sum_spec =~ s/(\d+)//) {
    push @cols_to_sum, $1;
    push @signs, '+';
  } else {
    print STDERR "Expected something like '4+8-5', got $col_sum_spec. Skipping this histogram.\n";
    return undef;
  }
  while ($col_sum_spec =~ s/([+-])(\d+)//) {
    push @cols_to_sum, $2;
    push @signs, $1;
  }

  # read lines store specified value if defined
  my $integer_data = 1;
  while (my $line = <$fh_in>) {
    next if($line =~ /^\s*#/);	# skip comments
    $line =~ s/^\s+//;		# remove initial whitespace
    my @columns = split(/\s+/, $line);

    # check whether the values in the specified columns are all defined and numbers
    my $all_look_like_numbers = 1;
    for my $cth (@cols_to_sum) {
      my $data_item = $columns[ $cth - 1 ] // undef;
      if (defined $data_item  and  looks_like_number($data_item)) {
	$integer_data = 0 if(! ($data_item =~ /^[+-]?\d+\z/) );
      } else {
	$all_look_like_numbers = 0; # specified quantity is not defined for this line
	last;
      }
    }

    # 
    if ($all_look_like_numbers == 1) {
      my $value_to_histogram = 0;
      while (my ($i, $column) = each @cols_to_sum ) {
	my $sign = $signs[$i];
	my $data_item = $columns[ $column - 1 ];
	if ($sign eq '+') {
	  $value_to_histogram += $data_item;
	} elsif ($sign eq '-') {
	  $value_to_histogram -= $data_item;
	} else {
	  print STDERR "sign should be + or -, but is: $sign . Skipping this histogram.\n";
	  return undef;
	}
      }
      $filecol_hdata->{$histogram_id}->add_value( $value_to_histogram );
      $filecol_hdata->{'pooled'}->add_value( $value_to_histogram );
    } else {
      $filecol_hdata->{$histogram_id}->add_value(undef);
      $filecol_hdata->{'pooled'}->add_value(undef);
    }
  }
  close($fh_in);
  return $integer_data;
}

sub read_and_union{
  my $data_file = shift;
  my $histogram_id = shift;
  my $col_spec = shift;
  my $filecol_hdata = shift;
  my @cols_to_histogram = split(/[uU]/, $col_spec);
  open my $fh_in, "<", $data_file or die "Couldn't open $data_file for reading.\n";
  my $integer_data = 1;
  while (my $line = <$fh_in>) {
    next if($line =~ /^\s*#/);	# skip comments
    $line =~ s/^\s+//;
    my @columns = split(/\s+/, $line);
    for my $cth (@cols_to_histogram ) {
      my $data_item = $columns[ $cth - 1 ];
      if (looks_like_number( $data_item ) ) {
	$integer_data = 0 if(! ($data_item =~ /^[+-]?\d+\z/) );
	$filecol_hdata->{$histogram_id}->add_value( $data_item );
	$filecol_hdata->{'pooled'}->add_value( $data_item );
      } else {
	$filecol_hdata->{$histogram_id}->add_value(undef);
	$filecol_hdata->{'pooled'}->add_value(undef);
      }
    }
  }
  close $fh_in;
  return $integer_data;
}

sub read_and_divide{ # divide the value in one column by the value in another
  my $data_file = shift;
  my $histogram_id = shift;
  my $col_spec = shift;
  my $filecol_hdata = shift;
  my @cols_to_use = split(/[\/]/, $col_spec);
  if (scalar @cols_to_use != 2) {
    print STDERR "$col_spec has problem. Only one divisor allowed. Skipping this histogram.\n";
    return undef;
  }
  my ($numer_col, $denom_col) = @cols_to_use;
  
  open my $fh_in, "<", $data_file or die "Couldn't open $data_file for reading.\n";
  while (my $line = <$fh_in>) {
    next if($line =~ /^\s*#/);	# skip comments
    $line =~ s/^\s+//;
    my @columns = split(/\s+/, $line);
    my ($numerator, $denominator) = ($columns[$numer_col-1] // undef, $columns[$denom_col-1] // undef);
    if (defined $numerator and defined $denominator and looks_like_number($numerator) and looks_like_number($denominator)) {
      my $value_to_histogram = ($denominator == 0)? undef : $numerator/$denominator;
    #  print "n, d, n/d:  $numerator $denominator  $value_to_histogram\n";
      $filecol_hdata->{$histogram_id}->add_value($value_to_histogram);
      $filecol_hdata->{'pooled'}->add_value($value_to_histogram);
    }else{
      $filecol_hdata->{$histogram_id}->add_value(undef);
      $filecol_hdata->{'pooled'}->add_value(undef);
    }
  }
  close $fh_in;
  return 0;
}

sub read_and_multiply{ # multiply the value in one column by the value in another
  my $data_file = shift;
  my $histogram_id = shift;
  my $col_spec = shift;
  my $filecol_hdata = shift;
  my @cols_to_use = split(/[\*]/, $col_spec);
  if (scalar @cols_to_use != 2) {
    print STDERR "$col_spec has problem. Only two factors allowed. Skipping this histogram.\n";
    return undef;
  }
  my ($factor1_col, $factor2_col) = @cols_to_use;

  open my $fh_in, "<", $data_file or die "Couldn't open $data_file for reading.\n";
  while (my $line = <$fh_in>) {
    next if($line =~ /^\s*#/);	# skip comments
    $line =~ s/^\s+//;
    my @columns = split(/\s+/, $line);
    my ($factor1, $factor2) = ($columns[$factor1_col-1] // undef, $columns[$factor2_col-1] // undef);
    if (defined $factor1 and defined $factor2 and looks_like_number($factor1) and looks_like_number($factor2)) {
      my $value_to_histogram = $factor1*$factor2;
    #  print "n, d, n*d:  $factor1 $factor2  $value_to_histogram\n";
      $filecol_hdata->{$histogram_id}->add_value($value_to_histogram);
      $filecol_hdata->{'pooled'}->add_value($value_to_histogram);
    }else{
      $filecol_hdata->{$histogram_id}->add_value(undef);
      $filecol_hdata->{'pooled'}->add_value(undef);
    }
  }
  close $fh_in;
  return 0;
}

sub read_and_log{	  # take the log of value in the column
  # if value is 0, should count as 'underflow'
  my $data_file = shift;
  my $histogram_id = shift;
  my $col_spec = shift;
  my $filecol_hdata = shift;
  my $col_to_use;
  if ($col_spec =~ /^log(\d+)/) {
    $col_to_use = $1;
  } else {
    print STDERR "$col_spec has problem. Should be something like log13 for log of value in col 13. Skipping this histogram.\n";
    return undef;
  }

  open my $fh_in, "<", $data_file or die "Couldn't open $data_file for reading.\n";
  while (my $line = <$fh_in>) {
    next if($line =~ /^\s*#/);	# skip comments
    $line =~ s/^\s+//;
    my @columns = split(/\s+/, $line);
    my $arg = ($columns[$col_to_use-1] // undef);
    #	print "col_spec:  $col_spec  arg:  $arg \n";
    if (defined $arg and looks_like_number($arg)) {
      if ($arg > 0) {
	my $value_to_histogram = log($arg)/log(10.0);
	#  print "n, d, n*d:  $factor1 $factor2  $value_to_histogram\n";
	$filecol_hdata->{$histogram_id}->add_value($value_to_histogram);
	$filecol_hdata->{'pooled'}->add_value($value_to_histogram);
      } elsif ($arg == 0) {
#	die;
	$filecol_hdata->{$histogram_id}->add_to_log0_count(1);
	$filecol_hdata->{'pooled'}->add_to_log0_count(1);
      } else {
	$filecol_hdata->{$histogram_id}->add_value(undef);
	$filecol_hdata->{'pooled'}->add_value(undef);
      }
    } else {
      $filecol_hdata->{$histogram_id}->add_value(undef);
      $filecol_hdata->{'pooled'}->add_value(undef);
    }
  }
  close $fh_in;
  return 0;
}

sub data_min_max{
  my $self = shift;
  return ($self->filecol_hdata->{'pooled'}->min(), $self->filecol_hdata->{'pooled'}->max());
}

__PACKAGE__->meta->make_immutable;

1;


# ###### unused ####################

# sub read_and_sum{
#   my $data_file = shift;
#   my $histogram_id = shift;
#   my $cols_to_sum = shift;
#   my $filecol_hdata = shift;
#   # my $hdata_obj = shift;
#   # my $pooled_hdata_obj = shift;
#   open my $fh_in, "<", $data_file or die "Couldn't open $data_file for reading.\n";
#   my $integer_data = 1;

#   #open my $fh_in, "<", $datafile or die "Couldn't open $datafile for reading.\n";
#   while (my $line = <$fh_in>) {
#     next if($line =~ /^\s*#/);	# skip comments
#     $line =~ s/^\s+//;
#     #   $line =~ s/\s+$//;
#     my @columns = split(/\s+/, $line);
#     # print "\n\nline: $line \n";
#     # print "cols to sum: ", join(", ", @$cols_to_sum), "\n";

#     my $all_look_like_numbers = 1;
#     for my $cth (@$cols_to_sum) {
#       my $data_item = $columns[ $cth - 1 ] // undef;
#       #  print "cth $cth  data_item $data_item \n";
      
#       if (defined $data_item  and  looks_like_number($data_item)) {
# 	$integer_data = 0 if(! ($data_item =~ /^[+-]?\d+\z/) );
#       } else {
# 	$all_look_like_numbers = 0;
# 	last;
#       }
#     }
#     #  print "all_loook_like_numbers: $all_look_like_numbers \n";
#     if ($all_look_like_numbers == 1) {
#       my $value_to_histogram = 0;
#       while (my ($i, $cth) = each @$cols_to_sum ) {

# 	my $data_item = $columns[ $cth - 1 ];
# 	#		print "cth $cth   data item:  $data_item \n";
# 	$value_to_histogram += $data_item;
#       }
#       #   print STDERR "histogram_id: [$histogram_id]\n";
#       $filecol_hdata->{$histogram_id}->
# 	# $hdata_obj->
# 	add_value( $value_to_histogram );
#       $filecol_hdata->{'pooled'}->
# 	# $pooled_hdata_obj->
# 	add_value( $value_to_histogram );
#     } else {
#       # $value_to_histogram = undef;
#     }
#   }
#   close($fh_in);
#   return $integer_data;
# }

# sub read_and_diff{
#   my $data_file = shift;
#   my $histogram_id = shift;
#   my @cols_to_diff = @{my $cols2diff = shift};
#   my $filecol_hdata = shift;
#   if (scalar @cols_to_diff != 2) {
#     print STDERR "Can only take the difference of 2 columns. Skipping histogram $histogram_id.\n";
#   }
#   print STDERR "In read_and_diff; cols to diff: ", join(", ", @cols_to_diff), "\n";
 
#   open my $fh_in, "<", $data_file or die "Couldn't open $data_file for reading.\n";
#   my $integer_data = 1;

#   while (my $line = <$fh_in>) {
#     next if($line =~ /^\s*#/);	       # skip comments
#     $line =~ s/^\s+//;		       # remove initial whitespace
#     my @columns = split(/\s+/, $line); # split on whitespace
#     my $data_item0 = $columns[$cols_to_diff[0] - 1];
#     my $data_item1 = $columns[$cols_to_diff[1] - 1];
#     if (looks_like_number($data_item0)  and  looks_like_number($data_item1)) {
#       if ($integer_data == 1) {
# 	$integer_data = 0 if(! ($data_item0 =~ /^[+-]?\d+\z/  and  $data_item1 =~ /^[+-]?\d+\z/) );
#       }
#       my $value_to_histogram = $data_item0 - $data_item1;
#       print STDERR "histogram_id: [$histogram_id]\n";
#       $filecol_hdata->{$histogram_id}->add_value( $value_to_histogram );
#       $filecol_hdata->{'pooled'}->add_value( $value_to_histogram );
#     } else {
#       #   $value_to_histogram = undef; # 
#     }
   
#   }
#   close($fh_in);
#   return $integer_data;
# }

