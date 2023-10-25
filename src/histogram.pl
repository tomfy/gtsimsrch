#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use List::Util qw (min max sum);

# use Math::GSL::SF  qw( :all );

my $y_plot_factor = 1.08;
my $y_plot_factor_log = 1.5;
my $relative_frame_thickness = 1.5; # the thickness of frame lines rel to histogram itself

use File::Basename 'dirname';
use Cwd 'abs_path';
my ( $bindir, $libdir );
BEGIN {     # this has to go in Begin block so happens at compile time
  $bindir =
    dirname( abs_path(__FILE__) ) ; # the directory containing this script
  $libdir = $bindir . '/../lib';
  $libdir = abs_path($libdir);	# collapses the bin/../lib to just lib
}
use lib $libdir;

use Histograms;
use Plot_params;

{				      # main
  unlink glob ".gnuplot*.stderr.log"; # to avoid accumulation of log files.
  my $lo_limit = 'auto';	      # 0;
  my $hi_limit = undef;
  my $binwidth = undef;
  my $persist = 0;
  my $do_plot = 1;
  my $log_y = 0;
  my $key_horiz_position = 'middle'; # 'right';
#  print "#### $key_horiz_position \n"; sleep(1);
  my $key_vert_position = 'top';
  my $data = undef;
  my $gnuplot_command = undef;
  my $line_width = 2;
  my $terminal = 'x11'; # (for gnuplot case) qt also works. Others?
  my $ymin = 0;
  my $ymax = "*";
  my $ymax_log = "*";
  my $ymin_log = 0.8;
  my $output_filename = 'histogram'; # default is to just send to screen.
  my $show_on_screen = 1;
  my $write_to_png = 0;
  my $interactive = undef;
  my $enhanced = 0;
  my $vline_position = undef;
  my $plot_title = undef;
  my $x_axis_label = undef;
  my $y_axis_label = undef;
  my $tight = 1;
  my $plot_width = 640; # pixels
  my $plot_height = 480; # pixels
  my $histogram_color = undef;
  # other options for plot range: 'strict', 'loose_positive', 'loose'
  my $graphics = 'gnuplot';	# alternative is 'gd' 

  GetOptions(
	     'data_input|input=s' => \$data,
	     'output_filename=s' => \$output_filename, # to send output to a (png) file, specify a filename

	     # control of binning, x and y ranges:
	     'low_limit|xmin=f' => \$lo_limit,
	     'hi_limit|xmax=f' => \$hi_limit,
	     'bw|binwidth|width=f' => \$binwidth,
	     'logy!' => \$log_y,
	     'ymax=f' => \$ymax,
	     'log_ymax=f' => \$ymax_log,
	     'tight!' => \$tight,
	     # how to plot (gnuplot or GD, and plot parameters)
	     'graphics=s' => \$graphics,
	     # whether to plot, and where (screen or file)
	     'plot!' => \$do_plot, # -noplot to suppress plot - just see histogram as text.
	     'screen!' => \$show_on_screen,
	     'png!' => \$write_to_png,
	     'interactive!' => \$interactive, # if true, plot and wait for further commands, else plot and exit

	     'width=f' => \$plot_width,
	     'height=f' => \$plot_height,

	     'linewidth|line_width|lw=f' => \$line_width, # line thickness for histogram
	     'color=s' => \$histogram_color,

	     'title=s' => \$plot_title,
	     'x_axis_label|x_label|xlabel=s' => \$x_axis_label,
	     'y_axis_label|y_label|ylabel=s' => \$y_axis_label,
	     'h_key|key_horiz_position=s' => \$key_horiz_position,
	     'v_key|key_vert_positioqn=s' => \$key_vert_position,

	     'vline_position=f' => \$vline_position,

	     # relevant to gnuplot
	     'terminal=s' => \$terminal, # x11, qt, etc.
	     'command=s' => \$gnuplot_command,
	     'enhanced!' => \$enhanced,
	    );

  if (lc $graphics eq 'gd') {
    print STDERR "GD graphics; setting terminal to only supported options: png \n" if(lc $terminal ne 'png');
    $terminal = 'png';
  } elsif (lc $graphics eq 'gnuplot') {
  } else {
    die "Graphics option $graphics is unknown. Allowed options are 'gnuplot' and 'gd'.\n";
  }

  print "vline position: ", $vline_position // "undef", "\n";
#  print "# key horiz position: $key_horiz_position \n"; sleep(2);
  # print "color: ", $histogram_color // "undef", "\n";
  # print "title: ", $plot_title // "undef", "\n";
  # print "x_axis_label: ", $x_axis_label // "undef", "\n";
  # print "y_axis_label: ", $y_axis_label // "undef", "\n";
  print "binwidth: ", $binwidth // "undef", "\n";

  $lo_limit = undef if($lo_limit eq 'auto'); # now default is 0.

  if (!defined $interactive) {
    $interactive = ($write_to_png)? 0 : 1;
  }

  $enhanced = ($enhanced)? 'enhanced' : 'noenhanced';

  print "files&columns to histogram: [$data] \n";

  my $histogram_obj = Histograms->new({
				       data_fcol => $data,
				       lo_limit => $lo_limit,
				       hi_limit => $hi_limit,
				       binwidth => $binwidth,
				       tight => $tight,
				      });
  $histogram_obj->bin_data();
  print "Max bin y: ", $histogram_obj->max_bin_y(), "\n";
  my $histogram_as_string = $histogram_obj->as_string();
  print "$histogram_as_string \n";

  print "graphics will use: ", lc $graphics, "\n";

  $ymax_log = $y_plot_factor_log*$histogram_obj->max_bin_y();
  $ymax = $y_plot_factor*$histogram_obj->max_bin_y();
  my $plot_params = Plot_params->new(
				     output_filename => $output_filename,
				     terminal => $terminal,
				     width => $plot_width,
				     height => $plot_height,
				     line_width => $line_width,
				     histogram_color => $histogram_color,
				     plot_title => $plot_title,
				     x_axis_label => $x_axis_label,
				     y_axis_label => $y_axis_label,
				     key_horiz_position => $key_horiz_position,
				     key_vert_position => $key_vert_position,

				     log_y => $log_y,
				     ymin => $ymin, ymax => $ymax,
				     ymin_log => $ymin_log, ymax_log => $ymax_log,

				     vline_position => $vline_position,
				    );
  #print "# line width: ", $plot_params->line_width(), "\n"; # sleep(4);

  if (lc $graphics eq 'gnuplot') {	# Use gnuplot

    
    use Graphics::GnuplotIF qw(GnuplotIF);
    my $plot_gnuplot = Graphics::GnuplotIF->new( persist => $persist, style => 'histeps');
    $plot_gnuplot->gnuplot_cmd("set terminal $terminal noenhanced linewidth $line_width size $plot_width, $plot_height");
    $plot_gnuplot->gnuplot_set_xlabel($x_axis_label) if(defined $x_axis_label);
    $plot_gnuplot->gnuplot_set_ylabel($y_axis_label) if(defined $y_axis_label);

    if ($log_y) {
      $plot_gnuplot->gnuplot_cmd('set log y');
      $plot_gnuplot->gnuplot_set_yrange($ymin_log, (defined $ymax_log)? $ymax_log : '*');
      set_arrow($plot_gnuplot, $vline_position, $ymin_log, $ymax_log) if(defined $vline_position);
    } else {
      if (defined $vline_position) {
	set_arrow($plot_gnuplot, $vline_position, $ymin, $ymax);
      }
      $plot_gnuplot->gnuplot_set_yrange(0, (defined $ymax)? $ymax : '*');
    }

    my $key_pos_cmd = 'set key ' . "$key_horiz_position  $key_vert_position";
    $plot_gnuplot->gnuplot_cmd($key_pos_cmd);
    $plot_gnuplot->gnuplot_cmd("set border lw $relative_frame_thickness"); # apparently this width is relative to that for the histogram lines.
    $plot_gnuplot->gnuplot_cmd('set mxtics');
    $plot_gnuplot->gnuplot_cmd('set tics out');
    $plot_gnuplot->gnuplot_cmd('set tics scale 2,1');
    $plot_gnuplot->gnuplot_cmd($gnuplot_command) if(defined $gnuplot_command);

    if ($write_to_png) {
      $plot_gnuplot->gnuplot_hardcopy($output_filename, " png $enhanced linewidth $line_width");
      $plot_gnuplot->gnuplot_cmd("set out $output_filename");
      plot_the_plot_gnuplot($histogram_obj, $plot_gnuplot);
      $plot_gnuplot->gnuplot_restore_terminal();
    }
    print "[$show_on_screen] [$terminal] [$do_plot]\n";
    if ($show_on_screen) {
      plot_the_plot_gnuplot($histogram_obj, $plot_gnuplot) if($do_plot);
    }
    if ($interactive) {
      #####  modify plot in response to keyboard commands: #####
      while (1) {		# loop to handle interactive commands.
	my $commands_string = <STDIN>; # command and optionally a parameter, e.g. 'x:0.8'
	last if(handle_interactive_command($histogram_obj, $plot_params,  $plot_gnuplot, $commands_string));
      }
    }
  }elsif (lc $graphics eq 'gd') { # Use GD
    use GD;
      plot_the_plot_gd($histogram_obj,
		       # {output_filename => $output_filename,
		       # 	width => $plot_width,
		       # 	height => $plot_height,
		       # 	linewidth => $line_width,
		       # 	histogram_color => $histogram_color,
		       # 	plot_title => $plot_title,
		       # 	x_axis_label => $x_axis_label,
		       # 	y_axis_label => $y_axis_label,
		       # 	key_horiz_position => $key_horiz_position,
		       # 	key_vert_position => $key_vert_position,

		       # 	vline_position => $vline_position}
		       $plot_params
		      );
  # my $vline_position = $parameters->{vline_position} // undef;
  # my $histogram_color = $parameters->{histogram_color} // undef; # default is black (then blue, green, red, ... if multiple histograms)
  # my $plot_title = $parameters->{plot_title} // undef;
  # my $x_axis_label = $parameters->{x_axis_label} // undef;
  # my $y_axis_label = $parameters->{y_axis_label} // undef;

      
    }
  else{
    die "Graphics option $graphics is unknown. Accepted options are 'gnuplot' and 'gd'\n";
  }
} # end of main
###########


sub plot_the_plot_gnuplot{
  my $histogram_obj = shift;
  my $plot_gnuplot_obj = shift;

  $plot_gnuplot_obj->gnuplot_set_xrange($histogram_obj->lo_limit(), $histogram_obj->hi_limit());
  my $bin_centers = $histogram_obj->filecol_hdata()->{pooled}->bin_centers();
  my @plot_titles = ();

  my @histo_bin_counts = ();
  while (my ($i_histogram, $plt) = each @{$histogram_obj->histograms_to_plot()} ) {
    if ($plt == 1) {
      print STDERR "adding histogram w index $i_histogram.\n";
      my $fcspec = $histogram_obj->filecol_specifiers()->[$i_histogram];
      my $hdata_obj = $histogram_obj->filecol_hdata()->{$fcspec};
      my $bincounts = $hdata_obj->bin_counts();
      push @plot_titles, $hdata_obj->label();
      push @histo_bin_counts, $bincounts; # push an array ref holding the bin counts ...
    }
  }

  $plot_gnuplot_obj->gnuplot_set_plot_titles(@plot_titles);
  $plot_gnuplot_obj->gnuplot_plot_xy($bin_centers, @histo_bin_counts);
}

sub set_arrow{
  my $the_plot = shift;
  my $x_pos = shift;
  my $y_min = shift;
  my $y_max = shift;
  print "arrow top: $y_max\n";
  $the_plot->gnuplot_cmd("unset arrow");
  $the_plot->gnuplot_cmd("set arrow nohead from $x_pos,$y_min to $x_pos,$y_max lw 0.75 dt '-'");
}


### GD plotting:
sub plot_the_plot_gd{
  my $histogram_obj = shift;
  my $parameters = shift;
  #my $parameters = shift; #
  my $output_filename = $parameters->output_filename();
  my $width = $parameters->width();
  my $height = $parameters->height();
  my $line_width = $parameters->line_width(); # line width for histogram, frame etc. are relative to this.
  my $histogram_color = $parameters->histogram_color() // undef; # default is black (then blue, green, red, ... if multiple histograms)
  my $plot_title = $parameters->plot_title() // undef;
  my $x_axis_label = $parameters->x_axis_label() // undef;
  my $y_axis_label = $parameters->y_axis_label() // undef;
  my $key_horiz_position = $parameters->key_horiz_position();
  my $key_vert_position = $parameters->key_vert_position();
  my $ymin = $parameters->ymin();
  my $ymax = $parameters->ymax();

  my $vline_position = $parameters->vline_position() // undef;

  my $filecol_hdata = $histogram_obj->filecol_hdata();
  open my $fhout, ">", "$output_filename";

  # print "vline position: ", $vline_position // "undef", "\n";
  # print "color: ", $histogram_color // "undef", "\n";
  # print "title: ", $plot_title // "undef", "\n";
  # print "x_axis_label: ", $x_axis_label // "undef", "\n";
  # print "y_axis_label: ", $y_axis_label // "undef", "\n";
  # sleep(5);

  my $char_width = 8; # these are the dimensions of GD's
  my $char_height = 16; # gdLargeFont characters, according to documentation.
  my $image = GD::Image->new($width, $height);
  my %colors = (black => $image->colorAllocate(0, 0, 0),
		white => $image->colorAllocate(255, 255, 255),
		blue => $image->colorAllocate(50, 80, 255),
		green => $image->colorAllocate(20,130,20),
		red => $image->colorAllocate(150,20,20));
  my $black = $colors{black};
  my @histogram_colors = ('black', 'blue', 'green', 'red');
  $image->filledRectangle(0, 0, $width-1, $height-1, $colors{white});

  my $margin = 28;		# margin width (in pixels)
  my $tick_length_pix = 6;
  my $x_axis_label_space = 3*$tick_length_pix + 2.5*$char_height;
  my $y_axis_label_space = 3*$tick_length_pix + (0.5 + length int($histogram_obj->max_bin_y()))*$char_width + ((defined $y_axis_label)? 2*$char_height : 0);

  my $frame_L_pix = $margin + $y_axis_label_space;
  my $frame_R_pix = $width - $margin;
  my $frame_T_pix = $margin;
  my $frame_B_pix = $height - ($margin + $x_axis_label_space);

  my $h_label_x_fraction = 0.1; # default histogram label horiz position is 'left'
  if($key_horiz_position eq 'middle'  or  $key_horiz_position eq 'center'){
    $h_label_x_fraction = 0.35;
  }elsif($key_horiz_position eq 'right'){
    $h_label_x_fraction = 0.6;
  }
  my $label_x_pix = (1.0 - $h_label_x_fraction)*$frame_L_pix + $h_label_x_fraction*$frame_R_pix;
  
  my $frame_line_thickness = $relative_frame_thickness*$line_width;
  $image->setThickness($frame_line_thickness);
  my $frame = GD::Polygon->new();
  $frame->addPt($frame_L_pix, $frame_B_pix);
  $frame->addPt($frame_L_pix, $frame_T_pix);
  $frame->addPt($frame_R_pix, $frame_T_pix);
  $frame->addPt($frame_R_pix, $frame_B_pix);
  $image->openPolygon($frame, $black);

  my $xmin = $histogram_obj->lo_limit();
  my $xmax = $histogram_obj->hi_limit();
 # my $ymin = 0;
 # my $ymax = $histogram_obj->max_bin_y()*$y_plot_factor;
  # log scale not implemented, so don't need $ymin_log, $ymax_log
  my $bin_width = $histogram_obj->binwidth();

  # add tick marks.
  #    on x axis:
  my $tick_x = 0;
  my $tick_spacing_x = 10*$bin_width;
  for (my $i=0; 1; $i++) {
    next if($tick_x < $xmin);
    last if($tick_x > $xmax);
    my $xpix = pix_pos($tick_x, $xmin, $xmax, $frame_L_pix, $frame_R_pix);
    if ($i%5 == 0) {
      $image->line($xpix, $frame_B_pix, $xpix, $frame_B_pix + 2*$tick_length_pix, $black);
      $image->string(gdLargeFont, $xpix - 0.5*(length $tick_x)*$char_width, $frame_B_pix + 3*$tick_length_pix, $tick_x, $black);
      $image->line($xpix, $frame_T_pix, $xpix, $frame_T_pix - 2*$tick_length_pix, $black);
    } else {
      $image->line($xpix, $frame_B_pix, $xpix, $frame_B_pix + $tick_length_pix, $black);
      $image->line($xpix, $frame_T_pix, $xpix, $frame_T_pix - $tick_length_pix, $black);
    }
    $tick_x += $tick_spacing_x;
  }
  # add a label to the x axis:
  if (defined $x_axis_label) {
    my $label_length = length $x_axis_label;
    my $xpix = pix_pos(0.5*($xmin + $xmax), $xmin, $xmax, $frame_L_pix, $frame_R_pix);
    $image->string(gdLargeFont,
		   $xpix - 0.5*$label_length*$char_width,
		   $frame_B_pix + 3*$tick_length_pix + $char_height,
		   $x_axis_label, $black);
  }

  # tick marks on y axis
  my $max_bin_count = $histogram_obj->max_bin_y();
  my $tick_y = 0;
  my $tick_spacing_y = tick_spacing($max_bin_count);
  my $max_tick_chars = 0; # max number of chars in the tick_y numbers
  for (my $i=0; 1; $i++) {
    next if($tick_y < $ymin);
    last if($tick_y > $ymax);
    my $ypix = pix_pos($tick_y, $ymin, $ymax, $frame_B_pix, $frame_T_pix);
    if ($i%5 == 0) {
      $image->line($frame_L_pix, $ypix, $frame_L_pix - 2*$tick_length_pix, $ypix, $black);
      $image->string(gdLargeFont, $frame_L_pix - 3*$tick_length_pix - (length $tick_y)*$char_width, $ypix-0.5*$char_height, $tick_y, $black);
      $image->line($frame_R_pix, $ypix, $frame_R_pix + 2*$tick_length_pix, $ypix, $black);
      $max_tick_chars = max($max_tick_chars, length $tick_y);
    } else {
      $image->line($frame_L_pix, $ypix, $frame_L_pix - $tick_length_pix, $ypix, $black);
      $image->line($frame_R_pix, $ypix, $frame_R_pix + $tick_length_pix, $ypix, $black);
    }
    $tick_y += $tick_spacing_y;
  }

    # add a label to the y axis:
  if (defined $y_axis_label) {
    my $label_length = length $y_axis_label;
    my $ypix = pix_pos(0.5*($ymin + $ymax), $ymin, $ymax, $frame_B_pix, $frame_T_pix);
    $image->stringUp(gdLargeFont, $frame_L_pix - (3*$tick_length_pix + ($max_tick_chars)*$char_width + 1.5*$char_height),
		   $ypix + 0.5*$label_length*$char_width,
		   $y_axis_label, $black);
  }

  ##########################
  ##  draw the histogram  ##
  ##########################
  my $histogram_line_thickness = $line_width;
  $image->setThickness($histogram_line_thickness);
  my @ids = keys %$filecol_hdata;
  my $n_histograms = (scalar @ids); # - 1; # subtract 1 to exclude the pooled histogram
 
#  for (my $i = 0; $i < $n_histograms; $i++) {
  while (my ($histogram_index, $plt) = each @{$histogram_obj->histograms_to_plot()} ) {
    if ($plt == 1) {
      print STDERR "adding histogram w index $histogram_index.\n";
      my $fcspec = $histogram_obj->filecol_specifiers()->[$histogram_index];
      my $hdata_obj = $histogram_obj->filecol_hdata()->{$fcspec};
      my $bincounts = $hdata_obj->bin_counts();
      # my $id = $ids[$i];
      # my $v = $filecol_hdata->{$id};
      # print "histogram i, id: $i  $id\n";
      my $hline = GD::Polygon->new(); # this will be the line outlining the histogram shape.
      my $bincenters = $hdata_obj->bin_centers();
      #   my $counts = $v->bin_counts();
      my $xpix = $frame_L_pix;
      my $ypix = $frame_B_pix;
      $hline->addPt($xpix, $ypix);
      while (my($i, $bcx) = each @$bincenters) {
	next if($bcx < $xmin  or  $bcx > $xmax); # exclude underflow, overflow
	my $bincount = $bincounts->[$i];
	$ypix = pix_pos($bincount, $ymin, $ymax, $frame_B_pix, $frame_T_pix);
	$xpix = pix_pos($bcx-0.5*$bin_width, $xmin, $xmax, $frame_L_pix, $frame_R_pix);
	$hline->addPt($xpix, $ypix);
	$xpix = pix_pos($bcx+0.5*$bin_width, $xmin, $xmax, $frame_L_pix, $frame_R_pix);
	$hline->addPt($xpix, $ypix);
      }				# end loop over histogram bins.
      $xpix = $frame_R_pix;
      $ypix = $frame_B_pix;
      $hline->addPt($xpix, $ypix);
      my $color = ($histogram_index == 0  and  defined $histogram_color)? $histogram_color : $histogram_colors[$histogram_index % 4];
      # print "[$histogram_index]  [$histogram_color]  [$color]\n";
      $image->unclosedPolygon($hline, $colors{$color});
      # my $label = $hdata_obj->label();
      if (defined $hdata_obj->label()  and  length $hdata_obj->label() > 0) {

	my $label_y_pix = 0.94*$frame_T_pix + 0.06*$frame_B_pix + $histogram_index*$char_height;
	$image->line($label_x_pix - 25, $label_y_pix + 0.5*$char_height, $label_x_pix - 5, $label_y_pix + 0.5*$char_height, $colors{$color});
	$image->string(gdLargeFont, $label_x_pix, $label_y_pix, $hdata_obj->label(), $black);
      }
    }
  }				# end loop over histograms
  if (defined $vline_position) {
    $image->setThickness(2);
   # my $black = $colors{black};
    $image->setStyle(
		     $black, $black, $black, $black, $black, $black, $black, $black, $black, $black, $black, $black,
		     $black, $black, $black, $black, $black, $black, $black, $black, $black, $black, $black, $black,
		     gdTransparent, gdTransparent, gdTransparent, gdTransparent, gdTransparent, gdTransparent,
		     gdTransparent, gdTransparent, gdTransparent, gdTransparent, gdTransparent, gdTransparent,
		     gdTransparent, gdTransparent, gdTransparent, gdTransparent);
    my $xpix = pix_pos($vline_position, $xmin, $xmax, $frame_L_pix, $frame_R_pix);
    $image->line($xpix, $frame_B_pix, $xpix, $frame_T_pix, gdStyled);
  }
  my $label = 'plot title';
  if(defined $plot_title){
  my $title_x_pix = 0.5*$frame_L_pix + 0.5*$frame_R_pix - (length $label)*$char_width;
  my $title_y_pix = 0.97*$frame_T_pix + 0.03*$frame_B_pix;
  $image->string(gdLargeFont, $title_x_pix, $title_y_pix, $plot_title, $black);
}

  binmode $fhout;
  print $fhout $image->png;
  close $fhout;
}

sub pix_pos{
  my $x = shift;
  my $xmin = shift;
  my $xmax = shift;
  my $low_edge = shift;		# the pixels corresponding to 
  my $hi_edge = shift;
  return $low_edge + ($x-$xmin)/($xmax-$xmin) * ($hi_edge-$low_edge);
}

sub tick_spacing{
  # put approx. 20 tick marks
  my $max_data = shift;
  my @spacing_options = (1,2,4,5,10);
  my $int_log10_max_data = int( log($max_data)/log(10) );
  my $z = $max_data/(10**$int_log10_max_data); # should be in range 1 <= $z < 10
  for my $sopt (@spacing_options) {
    if ($sopt > $z) {
      my $ts = $sopt*(10**$int_log10_max_data)/20;
      return $ts;
    }
  }
  print STDERR "### $max_data  $int_log10_max_data \n";
}


sub handle_interactive_command{ # handle 1 line of interactive command, e.g. r:4 or xmax:0.2;ymax:2000q
  my $histogram_obj = shift;
  my $the_plot_params = shift; # instance of Plot_params
   # my $plot_params = shift; # Plot_params->new( log_y => $log_y, ymin => $ymin, ymax => $ymax, ymin_log => $ymin_log, ymax_log => $ymax_log,
   # 			    vline_position => $vline_position, line_width => $line_width);
  my $the_gnuplot = shift; # $plot_gnuplot instance of Graphics::GnuplotIF
  my $commands_string = shift;
#  my $plot = $the_plot_params->plot_obj();
  my $log_y = $the_plot_params->log_y();
  my $ymin = $the_plot_params->ymin();
  my $ymax = $the_plot_params->ymax();
  my $ymin_log = $the_plot_params->ymin_log();
  my $ymax_log = $the_plot_params->ymax_log();
  my $vline_position = $the_plot_params->vline_position(); # 
  my $line_width = $the_plot_params->line_width();

  $commands_string =~ s/\s+$//g; # delete whitespace
  return 1 if($commands_string eq 'q');
  my @cmds = split(';', $commands_string);
  my ($cmd, $param) = (undef, undef);

  for my $cmd_param (@cmds) {
    if ($cmd_param =~ /^([^:]+):(.+)/) {
      ($cmd, $param) = ($1, $2);
      $param =~ s/\s+//g if(! $cmd =~ /^\s*xlabel\s*$/);
      print STDERR "cmd: [$cmd]  param: [$param]\n";
    } elsif ($cmd_param =~ /^(\S+)/) {
      $cmd = $1;
    }
    $cmd =~ s/\s+//g;

    if (defined $cmd) {
      $commands_string =~ s/\s+//g if($cmd ne 'xlabel');
      if ($cmd eq 'g') {
	$the_gnuplot->gnuplot_cmd('set grid');
      } elsif ($cmd eq 'll') {
	if ($log_y) {
	  $the_plot_params->log_y(0);
	  $the_gnuplot->gnuplot_cmd('unset log');
	  set_arrow($the_gnuplot, $vline_position, $ymin, $ymax) if(defined $vline_position);
	  $the_gnuplot->gnuplot_set_yrange($ymin, $ymax);
	} else {
	  $the_plot_params->log_y(1);
	  $the_gnuplot->gnuplot_cmd('set log y');
	  print STDERR "ll max_bin_y: ", $histogram_obj->max_bin_y(), " y_plot_factor: $y_plot_factor \n";
	  print STDERR "ll $ymin $ymax\n";
	  print "ymin,ymax,yminlog,ymaxlog: $ymin $ymax $ymin_log $ymax_log\n";
	  set_arrow($the_gnuplot, $vline_position, $ymin_log, $ymax_log) if(defined $vline_position);
	  $the_gnuplot->gnuplot_set_yrange($ymin_log, $ymax_log);
	}
      } elsif ($cmd eq 'ymax') {
	if (!$log_y) {
	  $ymax = $param;
	  set_arrow($the_gnuplot, $vline_position, $ymin, $ymax) if(defined $vline_position);
	  print "ymin,ymax,yminlog,ymaxlog: $ymin $ymax $ymin_log $ymax_log\n";
	  $the_gnuplot->gnuplot_set_yrange($ymin, $ymax);
	 # $the_plot_params->ymax($ymax);
	} else {
	  $ymax_log = $param;
	  set_arrow($the_gnuplot, $vline_position, $ymin_log, $ymax_log) if(defined $vline_position);
	  $the_gnuplot->gnuplot_set_yrange($ymin_log, $ymax_log);
	 # $the_plot_params->ymax_log($ymax_log);
	}
      } elsif ($cmd eq 'ymin') {
	if (!$log_y) {
	  $ymin = $param;
	  print "ymin,ymax,yminlog,ymaxlog: $ymin $ymax $ymin_log $ymax_log\n";
	  $the_gnuplot->gnuplot_set_yrange($ymin, $ymax);
	 # $the_plot_params->ymin($ymin);
	} else {
	  $ymin_log = $param;
	  $the_gnuplot->gnuplot_set_yrange($ymin_log, $ymax_log);
	 # $the_plot_params->ymin_log($ymin_log);
	}
      } elsif ($cmd eq 'x') {	# expand (or contract) x range.
	$histogram_obj->expand_range($param);
	$histogram_obj->bin_data();
      } elsif ($cmd eq 'bw') {	# set the bin width
	$histogram_obj->set_binwidth($param);
	$histogram_obj->bin_data();
      } elsif ($cmd eq 'lo' or $cmd eq 'low' or $cmd eq 'xmin') { # change low x-scale limit
	$histogram_obj->change_range($param, undef);
	$histogram_obj->bin_data();
      } elsif ($cmd eq 'hi' or $cmd eq 'xmax') { # change high x-scale limit
	$histogram_obj->change_range(undef, $param);
	$histogram_obj->bin_data();
      } elsif ($cmd eq 'c') {	# coarsen bins
	$histogram_obj->change_binwidth($param);
	$histogram_obj->bin_data();
	$ymax = $histogram_obj->max_bin_y()*$y_plot_factor;
	# $the_plot_params->ymax($ymax);
	$ymax_log = $histogram_obj->max_bin_y()*$y_plot_factor_log;
	# $the_plot_params->ymax_log($ymax_log);
	if ($log_y) {
	  $the_gnuplot->gnuplot_set_yrange($ymin_log, $ymax_log);
	  set_arrow($the_gnuplot, $vline_position, $ymin_log, $ymax_log) if(defined $vline_position);
	} else {
	  $the_gnuplot->gnuplot_set_yrange($ymin, $ymax);
	  set_arrow($the_gnuplot, $vline_position, $ymin, $ymax) if(defined $vline_position);
	}
      } elsif ($cmd eq 'r') {	# refine bins
	$histogram_obj->change_binwidth($param? -1*$param : -1);
	$histogram_obj->bin_data();
	$ymax = $histogram_obj->max_bin_y()*$y_plot_factor;
	# $the_plot_params->ymax($ymax);
	$ymax_log = $histogram_obj->max_bin_y()*$y_plot_factor_log;
	# $the_plot_params->ymax_log($ymax_log);
	if ($log_y) {
	  $the_gnuplot->gnuplot_set_yrange($ymin_log, $ymax_log);
	  set_arrow($the_gnuplot, $vline_position, $ymin_log, $ymax_log) if(defined $vline_position);
	} else {
	  $the_gnuplot->gnuplot_set_yrange($ymin, $ymax);
	  set_arrow($the_gnuplot, $vline_position, $ymin, $ymax) if(defined $vline_position);
	}
      } elsif ($cmd eq 'key') { # move the key (options are left, right, top, bottom)
	my $new_key_position = $param // 'left'; #
	$new_key_position =~ s/,/ /; # so can use e.g. left,bottom to move both horiz. vert. at once
	$the_gnuplot->gnuplot_cmd("set key $new_key_position");
      } elsif ($cmd eq 'xlabel') {
	$param =~ s/^\s+//;
	$param =~ s/\s+$//;
	$param =~ s/^([^'])/'$1/;
	$param =~ s/([^']\s*)$/$1'/;
	print STDERR "param: $param \n";
	$the_gnuplot->gnuplot_cmd("set xlabel $param");
      } elsif ($cmd eq 'export') {
	$param =~ s/'//g; # the name of the file to export to; the format will be png, and '.png' will be added to filename

	$the_gnuplot->gnuplot_hardcopy($param, " png linewidth $line_width");
	plot_the_plot_gnuplot($histogram_obj, $the_gnuplot);
	$the_gnuplot->gnuplot_restore_terminal();
      } elsif ($cmd eq 'off') {
	$histogram_obj->histograms_to_plot()->[$param-1] = 0;
      } elsif ($cmd eq 'on') {
	$histogram_obj->histograms_to_plot()->[$param-1] = 1;
      } elsif ($cmd eq 'cmd') {
	if ($param =~ /^\s*['](.+)[']\s*$/) { # remove surrounding single quotes if present
	  $param = $1;
	}
	$the_gnuplot->gnuplot_cmd("$param");
      }
      print STDERR "max_bin_y: ", $histogram_obj->max_bin_y(), " y_plot_factor: $y_plot_factor \n";
      plot_the_plot_gnuplot($histogram_obj, $the_gnuplot);

      $the_plot_params->ymin($ymin);
      $the_plot_params->ymin_log($ymin_log);
      $the_plot_params->ymax($ymax);
      $the_plot_params->ymax_log($ymax_log);
    }
  }
  return 0;
}
