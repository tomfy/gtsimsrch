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

has label => (		  # string describing quantity being clustered
	      isa => 'Str',
	      is => 'ro',
	      default => 'none specified'
	     );

has xs => (			# numbers to be put into 2 clusters
	   isa => 'ArrayRef[Num]',
	   is => 'rw',
	   required => 1,
	  );

has pow => ( # if this is a number y, cluster x^y, if it is 'log', cluster log(x)
	    isa => 'Str',
	    is => 'rw',
	    default => 'log',
	   );

has median_denom => (
		     isa => 'Num',
		     is => 'ro',
		     default => 3000,
		     # required => 1,
		    );

has txs => (			# transformed xs
	    isa => 'ArrayRef[Num]',
	    is => 'rw',
	    required => 0,
	   );

has minx => (
	     isa => 'Maybe[Num]',
	     is => 'rw',
	     default => 0.001, # to avoid -inf in case of log transformed values.
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
  my @txs = @xs;		# sort {$a <=> $b} @{$self->xs()};
  while ($txs[0] < 0) {		# shift the negatives away
    shift @txs;
  }
  # while ($txs[-1] > 1.1){ # pop away invalid data
  #   pop @txs;
  # }
  my $pow = $self->pow();
  if ($pow eq 'log') {
    my $minx = $self->minx // 1.0/($self->median_denom());
    print STDERR "minx $minx \n";
    #  exit;
    @txs = map(max($minx, $_), @xs);
    # my $small_limit = 1e-8;
    # my $xsmall = undef;
    # for my $x (@txs) { # find first (i.e. least) number >= $small_limit
    #   if ($x >= $small_limit) {
    # 	$xsmall = $x;
    # 	last;
    #   }
    # }
    # for my $x (@txs) {		# numbers < $xsmall get set to $xsmall
    #   if ($x < $xsmall) {
    # 	$x = $xsmall
    #   } else {
    # 	last;
    #   }
    # }
  }
  #  print STDERR "#  size of txs array: ", scalar @txs, "\n";

  if ($pow eq 'log') {
    @txs = map(log($_), @txs);
  } else {
    @txs = map($_**$pow, @txs);
  }
  $self->txs(\@txs);
}


sub one_d_2cluster{		# cluster 1dim data into 2 clusters
  my $self = shift;
  my $pow = $self->pow(); # cluster x**$pow ( or cluster log(x) if $pow eq 'log' )

  my $n_pts = scalar @{$self->txs()};
  # print STDERR "label: ", $self->label(), "  npts: $n_pts   pow: $pow  ";
  my ($km_n_L, $km_h_opt, $km_mom, $q) = $self->kmeans_2cluster();
  # print STDERR "$km_n_L  $km_h_opt $km_mom  $q  \n";
  my ($kde_n_L, $kde_h_opt, $min_kde_est, $kde_q) = $self->kde_2cluster($km_n_L-1);
  #  $km_h_opt = $km_mom; # maybe mean of means is better?
  if ($pow eq 'log') {
    $km_h_opt = exp($km_h_opt);
    $kde_h_opt = exp($kde_h_opt);
  } else {
    # print STDERR "pow, etc: $pow $km_h_opt   ";
    $km_h_opt = $km_h_opt**(1/$pow);
    $kde_h_opt = $kde_h_opt**(1/$pow);
    # print STDERR " $km_h_opt \n";
  }
  return ($n_pts, $km_n_L, $n_pts-$km_n_L, $km_h_opt, $q, 
	  $kde_n_L, $n_pts-$kde_n_L, $kde_h_opt, $kde_q);
}


sub kmeans_2cluster{		# divide into 2 clusters by minimizing
  # n_left*var_left + n_right*var_right
  # for N pts, just consider all N-1 possible ways of partitioning
  # into non-empty L and R sets with every value in L set < every value in R set.
  my $self = shift;
  my $txs = $self->txs();	# array ref of transformed values.
  my @xsqrs = map($_*$_, @$txs);
  # while (my ($i, $xx) = each @$txs){
  #     print STDERR "$i $xx  ", $xsqrs[$i], "\n";
  # 	 }
  my $h_opt = -1;
  my ($n, $sumx, $sumxsqr) = (scalar @$txs, sum(@$txs), sum(@xsqrs)); # sum(map($_*$_, @xs)));
  my ($n_left, $sumx_left, $sumxsqr_left) = (0, 0, 0);
  my ($n_right, $sumx_right, $sumxsqr_right) = ($n, $sumx, $sumxsqr);
  # print STDERR "ABC: $n_left  $sumx_left    $n_right  $sumx_right  $sumxsqr_right  \n";
  my $mean_of_means_opt = -1;
  my $v_opt = $sumxsqr_right - $sumx_right*$sumx_right/$n_right;
  my $n_left_opt = -1;
  # print STDERR "$v_opt $n_left_opt $mean_of_means_opt $h_opt\n";
  for my $x (@$txs[0 .. $#$txs-1]) {
    $n_left++; $n_right--;
    $sumx_left += $x; $sumx_right -= $x;
 
    $sumxsqr_left += $x*$x; $sumxsqr_right -= $x*$x;
    #   print STDERR "$sumx_left ", $sumx_left/$n_left, "  $sumxsqr_left       $sumx_right  ", ($sumx_right/$n_right)**2, "   $sumxsqr_right   ", $sumxsqr_right/$n_right, "  ", $x*$x, "\n";

    my $v = ($sumxsqr_left - $sumx_left*$sumx_left/$n_left) +  ($sumxsqr_right - $sumx_right*$sumx_right/$n_right);
    #  print STDERR "$x  $v  $v_opt $n_left_opt $mean_of_means_opt $h_opt\n";
    if ($v < $v_opt) {	    # if best so far record new best solution.
      $v_opt = $v;
      $n_left_opt = $n_left;
      $mean_of_means_opt =  0.5*($sumx_left/$n_left + $sumx_right/$n_right);
      $h_opt = 0.5*($x + $txs->[$n_left]);
      # print STDERR "$v_opt $n_left_opt $mean_of_means_opt $h_opt\n";
    }
    # my ($Lmean, $Rmean) = ($sumx_left/$n_left, $sumx_right/$n_right);
    # $mean_of_means = 0.5*($Lmean + $Rmean);
    #  print STDERR "$x  ",  $txs->[$n_left-1], "  ", $txs->[$n_left], "  $n_left  $sumx_left  $Lmean   $n_right  $sumx_right  $Rmean    $mean_of_means;\n";
    # if ($mean_of_means < $txs->[$n_left]  and  $mean_of_means >= $x) { # this is the place
    #   $h_opt = 0.5*($x + $txs->[$n_left]);
    # print STDERR " Left: $sumx_left  $n_left $Lmean  Right: $sumx_right  $n_right  $Rmean  \n";
    #   last;q
    # }
  }
  
  #  if(1){
  # my ($mean_left, $mean_right, $mean) = ($sumx_left/$n_left, $sumx_right/$n_right, $sumx/$n);
  # my $var_left = ($sumxsqr_left/$n_left - $mean_left**2);
  # my $var_right = ($sumxsqr_right/$n_right - $mean_right**2);
  # my $var = ($sumxsqr/$n - $mean**2);
  # my $q = sqrt($var_left + $var_right)/($mean_right - $mean_left);
  # # my $q1 = (($L90 - $Lmedian) + ($Rmedian - $R90))/($Rmedian - $Lmedian);
  #my $q = qqq($txs, $n_left);
  #  print STDERR "# $var_left $var_right $var   $q  $q1\n";
  #  getchar();
  #  }
  # print STDERR "XXX:  $n_left_opt $h_opt $mean_of_means_opt \n";
  return ($n_left_opt, $h_opt, $mean_of_means_opt, qqq($self->xs(), $n_left_opt, 0.05));
}

sub qqq{ # intended to be a measure of how well-separated the 2 clusters are.
  # essentially the ratio of the density in the L cluster peak to
  # the density in the region between the clusters
  # small value indicates good separation.
  
  my $xs = shift;		# all data pts.
  my $nL = shift;		# number in L cluster
  my $qile = shift // 0.05;
  my $n = scalar @$xs;

  my $nR = $n - $nL;		  # number in R cluster
  my $L25 = $xs->[int(0.25*$nL)]; # 25%ile of L cluster.
  my $L75 = $xs->[int(0.75*$nL)]; # 75%ile of L cluster.
  my $nLq = int($qile*$nL);
  my $L1mq = $xs->[$nL - $nLq]; # 1-$qile of L cluster is to left of this. i.e. $nLq pts of L cluster are to R of this.
  my $Rx = $xs->[$nL + $nLq]; # $nLq pts. of R cluster are to L of this
  # L peak density (of middle 50% of L cluster): 0.5*$nL / ($L75 - $L25)
  # valley density (region containing number of pts equal to 2*$qile fraction of L cluster
  # 2*$qile*$nL / ($Rx - $L1mq)
  my $density_ratio_valley_to_peak = 4*$qile * ($L75-$L25)/($Rx-$L1mq);
  return $density_ratio_valley_to_peak; # 1 - ($R10 - $L90)/($Rmedian - $Lmedian);
}

sub kde_2cluster{
  # now refine using kernel density estimation.
  my $self = shift;
  my $i_opt = shift; # look for min of kde in neighborhood of $xar->[$i_opt]
  my $kernel_width = shift // $self->kernel_width();
  my $txs = $self->txs();  # shift;
  # my $kernel_width = $self->kernel_width();
  my $n = scalar @$txs;
  my $n_left = $i_opt+1;
  my $n_right = $n - $n_left;
  if (!defined $kernel_width) { # consider 30% of pts near initial guess $i_opt
    $kernel_width = $self->n_pts_width($txs, $i_opt) * $self->width_factor(); # sqrt(2.0);
    $self->kernel_width($kernel_width);
  }
  
  # print STDERR "# in kde_2cluster. kernel width: $kernel_width \n";

  my $kde_i_opt = $i_opt;
  my $kde_x_est = $txs->[$i_opt];
  my $min_kde_est = kde($txs, $kde_x_est, $kernel_width, $i_opt);

  for (my $j = $i_opt; $j >= max(0, $i_opt-int($n_left/4)); $j--) { # starting at kmeans opt, search toward left for best kde.
    ($kde_i_opt, $kde_x_est, $min_kde_est, my $done) = $self->few_kdes($j, 3, $kde_i_opt, $kde_x_est, $min_kde_est);
    last if($done);
  }

  for (my $j = $i_opt+1; $j <= min($i_opt+int($n_right/4), scalar @$txs - 1); $j++) { # starting at kmeans opt, search toward right for best kde.
    ($kde_i_opt, $kde_x_est, $min_kde_est, my $done) = $self->few_kdes($j, 3, $kde_i_opt, $kde_x_est, $min_kde_est);
    last if($done);
  }

  my $kde_n_left = $kde_i_opt + 1;
  # my $q = qqq($txs, $kde_n_left);
  # print STDERR "kde $q \n";
  return ($kde_n_left, $kde_x_est, $min_kde_est, qqq($self->xs(), $kde_n_left, 0.05));
}


sub few_kdes{
  my $self = shift;

  my $i = shift;
  my $n_between = shift // 1;

  my $kde_i_opt = shift;
  my $kde_x_est = shift;
  my $min_kde_est = shift;

  my $txs = $self->txs();
  my $kernel_width = $self->kernel_width();

  my $done = 0;
  for my $k (0 .. $n_between-1) {
    my $eps = (0.5 + $k)/$n_between;
    my $x = (1.0 - $eps)*$txs->[$i] + $eps*$txs->[$i+1];
    my $kde_est = kde($txs, $x, $kernel_width, $i);
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
  # print STDERR "sufficient width: $sufficient_width \n";
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
  my $w = shift; # really the half-width, i.e. goes to zero at |x| = +- $width
  my $x = shift;
  # return (abs($x) >= $w)? 0 : (cos(PI*$x/$w) + 1);
  my $xx = PI*$x/(2.0*$w);
  return (abs($xx) < 3.0)? (1.0/(1.0 + $xx*$xx) - 0.1)/0.9 : 0.0;
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

# cluster by defining a quality of separation quantity and maximize
# or rather find range over which is quantity is > half max, define
# cluster boundary as some fraction through this
sub find_cluster_at_left{
  my $self = shift;
  my $f = shift // 0.5;
  my @txs = @{$self->txs()};	# array ref of transformed values.
  my $Qmin = 2;
  
  my $n_data_pts = scalar @txs;

  my $maxQ = -1;
  my $opt_nL = -1;
  my $optH = -1;
  my @nlhqs = ();
  my $start_peak = -1; # when Q goes above 1.5*$Qmin this gets set to $nL
  my $end_peak = -1; # when $start_peak > 0 and Q returns to below $Qmin this gets set to $nL
    my ($Ledge, $Redge) = (1, -1);
  for my $i (8..int($n_data_pts/4 -1)) { # loop over different values of the L cluster (which has 4*$i pts)
    my $nL = 4*$i;		# size of L cluster
    last if($nL >= 0.85*$n_data_pts);
    #last if($nL >= 10000);
    my $Lq1 = 0.5*($txs[$i - 1] + $txs[$i]); # 1st quartile of L cluster
    my $Lq3 = 0.5*($txs[3*$i - 1] + $txs[3*$i]); # 3rd quartile of L cluster

    my $nV = $i;
    my ($L78ths, $Rx);
    if ($i % 2 == 0) {		# $i even
      $L78ths = 0.5*($txs[7*$i/2 - 1] + $txs[7*$i/2]);
      $Rx =  0.5*($txs[9*$i/2 - 1] + $txs[9*$i/2]);
    } else {			# $i odd
      $L78ths = $txs[(7*$i-1)/2];
      $Rx = $txs[(9*$i-1)/2];
    }
    my $H= 0.5*($txs[$nL-1] + $txs[$nL]); # defines R edge of L cluster.
    #my $delta = $L78ths - $H; # half-width of valley
    #   my $Rxx = 2*$H  - $L78ths; # another way of defining R edge of valley.
    #   my $iR = ir(\@txs, $nL, $Rxx); # find $iR, s.t. $Rxx is between $txs[$iR], $txs[$iR+1] 
    #my $Rx = 0.5*($txs[17*$i -1 ] + $txs[17*$i]); # other size of valley
    #  print "$Rx $Rxx \n";
    #   $Rxx = $Rx;
    
    my $dL = 2*$i/($Lq3 - $Lq1);
    my $dV = $i/($Rx - $L78ths);
    my $Q =  $dL/$dV; #  2*($Rx - $L78ths)/($Lq3 - $Lq1);    #    4*($Rx - $Lq94)/($Lq3 - $Lq1);
    if ($Q > $maxQ) { # new best $Q
      $maxQ = $Q;
      $opt_nL = $nL;
      $optH = $H;		# 0.5*($txs[$nL -1 ] + $txs[$nL]);
    }
    #  my $H = 0.5*($txs[$nL -1 ] + $txs[$nL]);
  #    print STDERR "$H  $nL  $dL  $dV   $Lq3 $Lq1  $Rx $L78ths  $Q\n";
    if ($Q >= $Qmin) {
      push @nlhqs, [$nL, $H, $Q];
      $start_peak = $nL if($start_peak > 0  and  $Q >= 1.5*$Qmin);
    } else {			# $Q < $Qmin
      if ($start_peak > 0) {
	if ($end_peak < 0) {
	  $end_peak = $nL;
	} else {
	  if ($nL > 1.5*$end_peak) {
	    print "breaking out of loop. $start_peak, $end_peak, $opt_nL, $nL\n";
	    last;
	  }
	}
      }
    }
  #  print "$nL $start_peak $end_peak $opt_nL $maxQ \n";
  } # end of loop over possible L cluster sizes.
  # print STDERR "# start_peak, end_peak: $start_peak, $end_peak \n";
#  print STDERR "$maxQ $opt_nL $optH  ", scalar @nlhqs, "\n";
#  print "maxQ: $maxQ \n";
my $H_mid_half_max = -1;
if (scalar @nlhqs > 0) {
  ($Ledge, $Redge) = (1, -1);
  my ($nLmin, $nLmax) = (undef, undef);
  my $nabove = 0;
  for my $anlhq (@nlhqs) {
    my ($nl, $h, $q) = @$anlhq;
    #    print join(", ", @$anlhq), "\n";
    if ($q >= 0.5*$maxQ) {
      $nabove++;
      if ($h < $Ledge) {
	$Ledge = $h;
	$nLmin = $nl;
      }
      if ($h > $Redge) {
	$Redge = $h;
	$nLmax = $nl;
      }
    }
  }
  #    print "$nabove  $Ledge $Redge \n";
  if ($nabove >= 1) {
    $H_mid_half_max = 0.5*($Ledge+$Redge);
  }		     # otherwise leave $H_mid_half_max as -1;
} else {	     # no pts with $Q > 2 leave $H_mid_half_max as -1;
}
  # return $optH value at which $maxQ occurs,
  # and $Ledge, $Redge, the Left and Right sides of the > half max range.
  my $halfN = 16;
  my $N = 2*$halfN;
  my $Qsum = 0;
  my ($Qsum_opt, $mid_opt) = (-1, undef);
  for my $i (0..$N){
$Qsum += $nlhqs[$i]->[2];
}
  for(my $i=1; $i+$N < scalar @nlhqs; $i++){
    $Qsum += $nlhqs[$i+$N]->[2] - $nlhqs[$i-1]->[2];
    my $mid = $i+$halfN;
   # print STDERR $mid, "  ", join("  ", @{$nlhqs[$mid]}), "  ", $Qsum/($N+1), "\n";
    if($Qsum > $Qsum_opt){
      $Qsum_opt = $Qsum;
      $mid_opt = $mid;
    }
  }
  my $Qavg_opt = $Qsum_opt/($N+1);
  # print STDERR "### $mid_opt  $Qsum_opt \n";
  return ($optH, $maxQ, $nlhqs[$mid_opt]->[1], $Qavg_opt, $Ledge, $Redge); # $H_mid_half_max);
}

  #########################################################

#   sub ir{
#     my $txs = shift;
#     my $nL = shift;
#     my $Rxx = shift;
#     my $ir;
#     for ($ir = $nL; $ir < scalar @$txs; $ir++) {
#       last if ($txs->[$ir] <= $Rxx  and  $txs->[$ir+1] > $Rxx);
#     }
#     return $ir;
#   }

# # cluster by defining a (inverse) quality of separation quantity and minimize
# sub two_cluster_x{
#   my $self = shift;
#   my @txs = @{$self->txs()};	# array ref of transformed values.
#   my $a = 0.125;		# 
#   my $n = scalar @txs;

#   my $maxQ = -1;
#   my $opt_nL = -1;
#   my $optH = -1;
#   #for my $nL (24..5000){ #$n-10){ # $nL is the L cluster size being tested in one pass through loop
#   for (my $nL = 24; $nL <= 5000; $nL += 3) {
#     my @Lxs = @txs[0..$nL-1];	# pts in the L cluster
#     my @Rxs = @txs[$nL..$#txs]; # pts in the R cluster
#     my $nR = scalar @Rxs;
#     my $Lq1 = quantile(\@Lxs, 0.25); # 1st quartile of L cluster
#     my $Lq3 = quantile(\@Lxs, 0.75); # 3rd quartile of L cluster
#     my $Lq95 = quantile(\@Lxs, 1-$a);
#     #  $Lq95 = min($Lq95, $Lxs[$#Lxs]);
#     my $aR = $a*$nL/$nR;	# fraction of R cluster in the valley
  
#     last if($a*$nL >= ($nR-1));
#     my $vRx = quantile(\@Rxs, $aR); # other side of valley
#     my $dL = 0.5*$nL/($Lq3 - $Lq1); # L cluster interquartile density
#     # print STDERR "### $aR  $vRx  $Lq95 \n";
#     next if($a*$nL == 0  or  ($vRx - $Lq95) == 0);
#     my $dV = 2*$a*$nL/($vRx - $Lq95); # valley density
#     my $Q = $dL/$dV;
#     my $H = 0.5*($txs[$nL -1 ] + $txs[$nL]);
#     if ($Q > $maxQ) {
#       $maxQ = $Q;
#       $opt_nL = $nL;
#       $optH = 0.5*($txs[$nL -1 ] + $txs[$nL]);
#     }
#     #  my $H = 0.5*($txs[$nL -1 ] + $txs[$nL]);
#     #print STDERR "$H  $nL  $dL  $dV   $Lq3 $Lq1  $vRx $Lq95  $Q\n";
#   }
#   return ($optH, $maxQ);
# }



# sub quantile{		       # e.g. for $q = 0.5, returns the median
#   my $ar = shift;
#   my $q = shift;		# 0<$q<1
#   my @xs = @$ar;
#   my $f = $q * scalar @xs;
#   my $i = int($f - 0.5);
#   # print STDERR "##  $q  $f  $i  ", scalar @xs, "\n";
#   my $xq;
#   if ($i >= scalar @xs ) {
#     $xq = $xs[$#xs];
#   } else {
#     $xq = $xs[$i] + ($xs[$i+1] - $xs[$i])*($f - ($i+0.5));
#   }
#   return $xq;
# }

1;
