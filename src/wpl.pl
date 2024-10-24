#!/usr/bin/perl -w
use strict;

# mean and mode of number of words per line (whitespace separated)
my $lines_then_words = shift // 1;
my $verbose = shift // 1;
my $n_words_to_print = shift // undef;

my %nwords_nlines = ();
my $nlines = 0;
my $total_nwords = 0;
while(<>){
  my @cols = split(" ", $_);
  my $nwords = scalar @cols;
  if(defined $n_words_to_print  and  $nwords == $n_words_to_print){
    print; # "$nwords [$_]\n";
  }
  $nwords_nlines{$nwords}++;
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

if ($verbose) {
  for my $nw (@sorted_nwords) {
    print "# there are ", $nwords_nlines{$nw}, " lines with $nw words.\n";
  }
}
