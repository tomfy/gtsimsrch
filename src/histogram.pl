#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use List::Util qw (min max sum);

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

{				      # main
  unlink glob ".gnuplot*.stderr.log"; # to avoid accumulation of log files.
  my $lo_limit = 'auto';	      # 0;
  my $hi_limit = undef;
  my $binwidth = undef;
  my $persist = 0;
  my $do_plot = 1;
  my $log_y = 0;
  my $key_horiz_position = 'center'; # 'right';
  #  print "#### $key_horiz_position \n"; sleep(1);
  my $key_vert_position = 'top';
  my $data = undef;
  my $gnuplot_command = undef;
  my $line_width = 2;
  my $terminal = 'x11';	   # (for gnuplot case) qt also works. Others?
  my $ymin = 0;
  my $ymax = "*";
  my $ymax_log = "*";
  my $ymin_log = 0.8;
  my $output_filename = 'histogram.png'; # default is to just send to screen.
  my $show_on_screen = 1;
  my $write_to_png = 0;
  my $interactive = undef;
  my $enhanced = 0;
  my $vline_position = undef;
  my $plot_title = undef;
  my $x_axis_label = undef;
  my $y_axis_label = undef;
  my $tight = 1;     # choose x range so as to include all data points
  my $plot_width = 640;		# pixels
  my $plot_height = 480;	# pixels
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
	     'v_key|key_vert_position=s' => \$key_vert_position,

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

  $lo_limit = undef if($lo_limit eq 'auto'); # now default is 0.

  if (!defined $interactive) {
    $interactive = ($write_to_png)? 0 : 1;
  }

  $enhanced = ($enhanced)? 'enhanced' : 'noenhanced';

  print "files&columns to histogram: [$data] \n";

  my $histograms_obj = Histograms->new({
					data_fcol => $data,
					lo_limit => $lo_limit,
					hi_limit => $hi_limit,
					binwidth => $binwidth,
					tight => $tight,
				       });
  $histograms_obj->bin_data();
  print "Max bin y: ", $histograms_obj->max_bin_y(), "\n";
  my $histogram_as_string = $histograms_obj->as_string();
  print "$histogram_as_string \n";

  print "graphics will use: ", lc $graphics, "\n";

  $ymax_log = $y_plot_factor_log*$histograms_obj->max_bin_y();
  $ymax = $y_plot_factor*$histograms_obj->max_bin_y();

  if (lc $graphics eq 'gnuplot') { # Use gnuplot
    use Gnuplot_plot;
    my $gnuplot_plot = Gnuplot_plot->new({
					  persist => $persist,
					  width => $plot_width, height => $plot_height,
					  xmin => $histograms_obj->lo_limit, xmax => $histograms_obj->hi_limit,
					  ymin => $ymin, ymax => $ymax,
					  ymin_log => $ymin_log, ymax_log => $ymax_log,
					  x_axis_label => $x_axis_label,
					  y_axis_label => $y_axis_label,
					  line_width => $line_width,
					  color => $histogram_color,
					  relative_frame_thickness => $relative_frame_thickness,
					  key_horiz_position => $key_horiz_position,
					  key_vert_position => $key_vert_position,
					  histograms => $histograms_obj,
					  vline_position => $vline_position,
					  output_filename => $output_filename,
					 });
    $gnuplot_plot->draw_histograms();

    if ($interactive) {
      #####  modify plot in response to keyboard commands: #####
      while (1) {		# loop to handle interactive commands.
	my $commands_string = <STDIN>; # command and optionally a parameter, e.g. 'x:0.8'
	my $done = $gnuplot_plot->handle_interactive_command($commands_string);
	#print "done with interactive commands? $done \n";
	last if($done);
      }
    }
    #}
  } elsif (lc $graphics eq 'gd') { # Use GD
    use GD_plot;
    my $gd_plot = GD_plot->new({
				persist => $persist,
				width => $plot_width, height => $plot_height,
				xmin => $histograms_obj->lo_limit, xmax => $histograms_obj->hi_limit,
				ymin => $ymin, ymax => $ymax,
				x_axis_label => $x_axis_label,
				y_axis_label => $y_axis_label,
				line_width => $line_width,
				color => $histogram_color,
				relative_frame_thickness => $relative_frame_thickness,
				key_horiz_position => $key_horiz_position,
				key_vert_position => $key_vert_position,
				histograms => $histograms_obj,
				vline_position => $vline_position,
				output_filename => $output_filename,
			       });
    $gd_plot->draw_histograms();
    $gd_plot->draw_vline($vline_position, 'black');
    print STDERR "unlink output_filename : ", $gd_plot->output_filename, "\n";
    unlink $gd_plot->output_filename;
    print STDERR "open output_filename : ", $gd_plot->output_filename, "\n";
    open my $fhout, ">", $gd_plot->output_filename;
    binmode $fhout;
    print $fhout $gd_plot->image->png;
    close $fhout;

    if ($interactive) {
      #####  modify plot in response to keyboard commands: #####
      while (1) {		# loop to handle interactive commands.
	my $commands_string = <STDIN>; # command and optionally a parameter, e.g. 'x:0.8'
	my $done = $gd_plot->handle_interactive_command($commands_string);
	last if($done);
      }
    }
    } else {
      die "Graphics option $graphics is unknown. Accepted options are 'gnuplot' and 'gd'\n";
    }
    print "Exiting histogram.pl\n";
  }
  ###########  end of main  ##############
