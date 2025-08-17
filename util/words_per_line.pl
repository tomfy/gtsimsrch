#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw(min max sum shuffle);
# mean and mode of number of words per line (whitespace separated)
my $lines_then_words = 1;
my $output_sample_line = 0;
my $n_words_to_print = undef;
my $n_lines_max = 1000000;
my $split_char = 'tab'; # set to anything else (with e.g. -split space ) to split on whitespace.
my $include_comments = 0;

GetOptions(
	   'nlines!' => \$lines_then_words, # 1 -> output in order of number of lines with each number of words
	   'sample_line!' => \$output_sample_line, # print the first line with each number of words
	   'words_out=i' => \$n_words_to_print,
	   'max_lines=i' => \$n_lines_max,
	   'split=s' => \$split_char,
	   'comments!' => \$include_comments,
	   );

my $skip_include_comments = ($include_comments)? "Including" : "Skipping";
print "# $skip_include_comments comments (i.e. lines with # as first non-whitespace character.)\n";
print "# Splitting on:  ", (($split_char eq 'tab')? 'tab' : 'whitespace'), "\n";

my %nwords_nlines = ();
my %nwords_firstline = ();
my $nlines = 0;
my $total_nwords = 0;
while(my $line = <>){
  next if(!$include_comments and $line =~ /^\s*#/);
  chomp($line);
  my @cols = ($split_char eq 'tab')? split("\t", $line) : split(" ", $line);
  # while(my ($i, $v) = each @cols){
  #   print "$i  [$v]\n";
  # }exit;
  my $nwords = scalar @cols;
  if(defined $n_words_to_print  and  $nwords == $n_words_to_print){
    print $line; # "$nwords [$line]\n";
  }
  $nwords_nlines{$nwords}++;
  if($nwords_nlines{$nwords} == 1){
    $nwords_firstline{$nwords} = $line;
  }
  $total_nwords += $nwords;
  $nlines++;
  last if($nlines >= $n_lines_max);
}

my $mean_nwords = ($nlines > 0)? $total_nwords/$nlines : '--';
printf("# lines: %ld ;  words: %ld ;  mean words per line: %8.2f\n", $nlines, $total_nwords, $mean_nwords);
my @sorted_nwords = ($lines_then_words)?
  sort {my $res = $nwords_nlines{$b} <=> $nwords_nlines{$a}; return ($res == 0)? ($b <=> $a)  : $res; } keys %nwords_nlines
  :
  sort {my $res = ($b <=> $a); return ($res == 0)? $nwords_nlines{$b} <=> $nwords_nlines{$a} : $res; } keys %nwords_nlines
  ;
my $str = ($lines_then_words)?
  "# most frequent number of words per line:  " . $sorted_nwords[0] . " ; which occurs on " . $nwords_nlines{$sorted_nwords[0]} . " lines."
  :
  "# max number of words per line:  " . $sorted_nwords[0] . " ; which occurs on " . $nwords_nlines{$sorted_nwords[0]} . " lines.";
print "$str\n";


for my $nw (@sorted_nwords) {
  print "# there are ", $nwords_nlines{$nw}, " lines with $nw words.\n";
  if ($output_sample_line){
    print "   ", $nwords_firstline{$nw}, "\n";
  }
}
 print "\n";
