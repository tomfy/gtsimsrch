package Cluster1d;
use strict;
use warnings;
use Moose;
use namespace::autoclean;
use Carp;
use Math::Trig;			# so can use pi() for pi.
use List::Util qw(min max sum);

use constant PI => pi();


# for clustering a set of non-negative numbers into 2 clusters.
# negative numbers are interpreted as invalid and excluded.

has label => ( # string describing quantity being clustered
	      isa => 'Str',
	      is => 'ro',
	      default => 'none specified'
	      );

has xs => ( # numbers to be put into 2 clusters
	   isa => 'ArrayRef[Num]',
	   is => 'rw',
	   required => 1,
	  );

has pow => ( # if this is a number y, cluster x^y, if it is 'log', cluster log(x)
	    isa => 'Str',
	    is => 'rw',
	    default => 'log',
	   );

has txs => (			# transformed xs
	    isa => 'ArrayRef[Num]',
	    is => 'rw',
	    required => 0,
	   );



has n_pts_in_kernel_width => (
			      isa => 'Int',
			      is => 'ro',
			      default => 16,
			     );

has width_factor => (
		     isa => 'Num',
		     is => 'ro',
		     default => 0.7, # sub { 1.0/sqrt(2.0) },
		    );

has kernel_width => (
		     isa => 'Maybe[Num]',
		     is => 'rw',
		     default => undef,
		    );




sub BUILD{ # for clustering values in range [0,1]; values outside are invalid - skip
  my $self = shift;
#  print STDERR "clustered quantity: ", $self->label(), "\n";
#  print STDERR join(" ", $self->xs()->[0..20]), "\n";
  my @xs = sort {$a <=> $b} @{$self->xs()};
  $self->xs(\@xs);
  my @txs = sort {$a <=> $b} @{$self->xs()};
  while ($txs[0] < 0) {		# shift the negatives away
    shift @txs;
  }
  while ($txs[-1] > 1.1){ # pop away invalid data
    pop @txs;
  }

  my $small_limit = 1e-8;
  my $xsmall = undef;
  for my $x (@txs) {  # find first (i.e. least) number >= $small_limit
    if ($x >= $small_limit) {
      $xsmall = $x;
       last;
    }
  }
  for my $x (@txs) {		# numbers < $xsmall get set to $xsmall
    if ($x < $xsmall) {
      $x = $xsmall
    } else {
      last;
    }
  }
#  print STDERR "#  size of txs array: ", scalar @txs, "\n";
  my $pow = $self->pow();
  if ($pow eq 'log') {
    @txs = map(log($_), @txs);
  } else {
    @txs = map($_**$pow, @txs);
  }
  $self->txs(\@txs);
}


sub one_d_2cluster{ # cluster 1dim data into 2 clusters
  my $self = shift;
  my $pow = $self->pow();	# cluster x**$pow (or log(x)

  my $n_pts = scalar @{$self->txs()};
  my ($km_n_L, $km_h_opt, $km_mom, $q) = $self->kmeans_2cluster();
  my ($kde_n_L, $kde_h_opt, $min_kde_est) = $self->kde_2cluster($km_n_L-1);
  if ($pow eq 'log') {
    $km_h_opt = exp($km_h_opt);
    $kde_h_opt = exp($kde_h_opt);
  } else {
   # print STDERR "pow, etc: $pow $km_h_opt   ";
    $km_h_opt = $km_h_opt**(1/$pow);
    $kde_h_opt = $kde_h_opt**(1/$pow);
  #  print STDERR " $km_h_opt \n";
  }
  return ($n_pts, $km_n_L, $n_pts-$km_n_L, $km_h_opt, $q, 
	  $kde_n_L, $n_pts-$kde_n_L, $kde_h_opt);
}


sub kmeans_2cluster{ # divide into 2 clusters by finding dividing value h s.t. h = 1/2 * (<x>_<h + <x>_>h)
  # however instead of guessing some initial cluster centers and iteratively refining, just 
  # consider break points with 1, 2, 3, etc. pts in the L-hand cluster,
  # until mean of means (i.e. mean of the mean of L cluster, and mean of R cluster)
  # lies between the two clusters.
  my $self = shift;
  my $xs = $self->txs(); # array ref of transformed values.
  my @xsqrs = map($_*$_, @$xs);
  my $h_opt;
  my ($n, $sumx, $sumxsqr) = (scalar @$xs, sum(@$xs), sum(@xsqrs)); # sum(map($_*$_, @xs)));
  my ($n_left, $sumx_left, $sumxsqr_left) = (0, 0, 0);
  my ($n_right, $sumx_right, $sumxsqr_right) = ($n, $sumx, $sumxsqr);

  my $mean_of_means;
  for my $x (@$xs[0 .. $#$xs-1]) {
    $n_left++; $n_right--;
    $sumx_left += $x; $sumx_right -= $x;
    $sumxsqr_left += $x*$x; $sumxsqr_right -= $x*$x;
 #   print STDERR "$sumx_left ", $sumx_left/$n_left, "  $sumxsqr_left       $sumx_right  ", ($sumx_right/$n_right)**2, "   $sumxsqr_right   ", $sumxsqr_right/$n_right, "  ", $x*$x, "\n";
    $mean_of_means = 0.5*($sumx_left/$n_left + $sumx_right/$n_right);
    if ($mean_of_means < $xs->[$n_left]  and  $mean_of_means >= $x) { # this is the place
      $h_opt = 0.5*($x + $xs->[$n_left]);
      last;
    }
  }
  if(1){

  # my ($mean_left, $mean_right, $mean) = ($sumx_left/$n_left, $sumx_right/$n_right, $sumx/$n);
  # my $var_left = ($sumxsqr_left/$n_left - $mean_left**2);
  # my $var_right = ($sumxsqr_right/$n_right - $mean_right**2);
  # my $var = ($sumxsqr/$n - $mean**2);
  # my $q = sqrt($var_left + $var_right)/($mean_right - $mean_left);
  # # my $q1 = (($L90 - $Lmedian) + ($Rmedian - $R90))/($Rmedian - $Lmedian);
  #my $q = qqq($xs, $n_left);
#  print STDERR "# $var_left $var_right $var   $q  $q1\n";
  #  getchar();
}
  return ($n_left, $h_opt, $mean_of_means, qqq($xs, $n_left));
}

sub qqq{
  my $xs = shift;
  my $nL = shift;
  my $n = scalar @$xs;
  my $nR = $n - $nL;
    my $Lmedian = $xs->[int(0.5*$nL)];
  my $Rmedian = $xs->[$n-1 - int(0.5*$nR)];
  my $L90 = $xs->[int(0.9*$nL)];
  my $R90 = $xs->[$n-1 - int(0.9*$nR)];
  return (($L90 - $Lmedian) + ($Rmedian - $R90))/($Rmedian - $Lmedian);
}

sub kde_2cluster{
  # now refine using kernel density estimation.
  my $self = shift;
  my $i_opt = shift; # look for min of kde in neighborhood of $xar->[$i_opt] 
 my $xs = $self->txs(); # shift;
  my $kernel_width = $self->kernel_width();
  my $n = scalar @$xs;
  my $n_left = $i_opt+1;
  my $n_right = $n - $n_left;
  if (!defined $kernel_width) { # consider 30% of pts near initial guess $i_opt
    $kernel_width = $self->n_pts_width($xs, $i_opt) * $self->width_factor(); # sqrt(2.0);
    $self->kernel_width($kernel_width);
  }
  # print STDERR "# in kde_2cluster. kernel width: $kernel_width \n";

  my $kde_i_opt = $i_opt;
  my $kde_x_est = $xs->[$i_opt];
  my $min_kde_est = kde($xs, $kde_x_est, $kernel_width, $i_opt);

  for (my $j = $i_opt; $j >= max(0, $i_opt-int($n_left/4)); $j--) { # starting at kmeans opt, search toward left for best kde.
    ($kde_i_opt, $kde_x_est, $min_kde_est, my $done) = $self->few_kdes($j, 3, $kde_i_opt, $kde_x_est, $min_kde_est);
    last if($done);
  }

  for (my $j = $i_opt+1; $j <= min($i_opt+int($n_right/4), scalar @$xs - 1); $j++) { # starting at kmeans opt, search toward right for best kde.
    ($kde_i_opt, $kde_x_est, $min_kde_est, my $done) = $self->few_kdes($j, 3, $kde_i_opt, $kde_x_est, $min_kde_est);
    last if($done);
  }

  my $kde_n_left = $kde_i_opt + 1;
  # my $q = qqq($xs, $kde_n_left);
  # print STDERR "kde $q \n";
  return ($kde_n_left, $kde_x_est, $min_kde_est);
}


sub few_kdes{
  my $self = shift;

  my $i = shift;
  my $n_between = shift // 1;

  my $kde_i_opt = shift;
  my $kde_x_est = shift;
  my $min_kde_est = shift;

  my $xs = $self->txs();
  my $kernel_width = $self->kernel_width();

  my $done = 0;
  for my $k (0 .. $n_between-1) {
    my $eps = (0.5 + $k)/$n_between;
    my $x = (1.0 - $eps)*$xs->[$i] + $eps*$xs->[$i+1];
    my $kde_est = kde($xs, $x, $kernel_width, $i);
    if ($kde_est < $min_kde_est) {
      $min_kde_est = $kde_est;
      $kde_x_est = $x;
      $kde_i_opt = $i;
    } elsif ($kde_est > 10*$min_kde_est) {
      $done = 1;		# last;
    }
    # print STDERR "$i  $x   $kde_est $min_kde_est\n";
  }
  return ($kde_i_opt, $kde_x_est, $min_kde_est, $done);
}


sub n_pts_width{ # in vicinity of $xs[$iguess], find interval width needed to guarantee including $n_points points
  my $self = shift;
  my $xs = shift;
  my $i_guess = shift; # initial guess - look in neighborhood of this index.
  my $n_points = $self->n_pts_in_kernel_width();

  my $n = scalar @$xs;
  my $eps = 0.3;
  my $istart = int( (1-$eps)*$i_guess + $eps*0 );
  my $iend = int( (1-$eps)*$i_guess + $eps*($n-1) );
  # print STDERR "# in n_pts_width n_points, istart, iend: $n_points  $istart $iend \n";
  my $sufficient_width = -1;
  for my $i ($istart .. $iend-$n_points) {
    my $width = $xs->[$i+$n_points-1] - $xs->[$i];
    if ($width > $sufficient_width) {
      $sufficient_width = $width;
    }
  }
  return $sufficient_width;
}


### non-methods:

sub kde{
  my $xs = shift;	     # array ref of sorted reals (low-to-hi)
  my $x = shift;	     # evaluate the kde at this x
  my $w = shift;	     # half-width of kernel
  my $i = shift // undef;    # largest index such that $xs->[$i] <= $x

  my $kde_sum = 0;

  for (my $j=$i; $j >= 0; $j--) { # go to smaller and smaller values of x until further than $w away
    my $arg = abs( $xs->[$j] - $x );
    last if ($arg >= $w);
    $kde_sum += kernel($w, $arg);
  }
  for (my $j=$i+1; $j < scalar @$xs; $j++) { # go to larger and larger values of x until further than $w away
    my $arg = abs( $xs->[$j] - $x );
    last if ($arg >= $w);
    $kde_sum += kernel($w, $arg);
  }
  return $kde_sum;
}

sub kernel{			# 2 at x=0, 0 at |x| >= w
  my $w = shift;
  my $x = shift;
  return (abs($x) >= $w)? 0 : (cos(PI*$x/$w) + 1)
}



sub jenks_2cluster{ # divide into 2 clusters using jenks natural breaks
  my $xs = shift;   # array ref of real data values
  # my @xs = sort {$a <=> $b} @$xar; # @xs is sorted small to large.

  my ($n_left_opt, $LR_max) = (1, -10000.0);
  my ($n, $sumx, $sumxsqr) = (scalar @$xs, sum(@$xs), 0); # sum(map($_*$_, @xs)));
  my ($n_left, $sumx_left, $sumxsqr_left) = (0, 0, 0);
  my ($n_right, $sumx_right, $sumxsqr_right) = ($n, $sumx, $sumxsqr);

  my $n_pts = scalar @$xs;
  $LR_max *= $n_pts;
  for my $x (@$xs[0 .. $#$xs-1]) {
    $n_left++; $n_right--;
    $sumx_left += $x; $sumx_right -= $x;
    #    $sumxsqr_left += $x*$x; $sumxsqr_right -= $x*$x;
    my $L = ($sumx_left)**2/$n_left;   # - $sumxsqr_left;
    my $R = ($sumx_right)**2/$n_right; # - $sumxsqr_right;
    my $LR = ($L + $R);
    if ($LR > $LR_max) {
      $LR_max = $LR;
      $n_left_opt = $n_left;
    }
  }
  my $h_opt = 0.5*($xs->[$n_left_opt-1] + $xs->[$n_left_opt]);
  return ($n_left_opt, $h_opt, $LR_max);
}


1;
