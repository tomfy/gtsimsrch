#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw(min max sum shuffle);
# mean and mode of number of words per line (whitespace separated)
my $lines_then_words = 1;
my $output_sample_line = 0;
my $n_words_to_print = undef;

GetOptions(
	   'nlines!' => \$lines_then_words, # 1 -> output in order of number of lines with each number of words
	   'sample_line!' => \$output_sample_line, # print the first line with each number of words
	   'words_out=i' => \$n_words_to_print,
	   );

my %nwords_nlines = ();
my %nwords_firstline = ();
my $nlines = 0;
my $total_nwords = 0;
while(my $line = <>){
  my @cols = split(" ", $line);
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
    print "   ", $nwords_firstline{$nw};
  }
}
