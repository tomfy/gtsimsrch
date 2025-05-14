#!/usr/bin/perl -w
use strict;

my $file1 = shift;
my $file2 = shift;
my $offset1 = shift // 0; # add to 3 to get column of first parent in file 1
my $offset2 = shift // 0; # add to 3 to column of first parent in file 2
my $p1column_file1 = 3 + $offset1;
$p1column_file1--;
my $p1column_file2 = 3 + $offset2;
$p1column_file2--;
my $valcolumnstr = shift // '7,8,9,10,11';
my @valcolumns_file1 = split(",", $valcolumnstr);
@valcolumns_file1 = map($_ + $offset1, @valcolumns_file1);
print "# columns: ", join("  ", @valcolumns_file1), "  of $file1 will be compared to \n";
@valcolumns_file1 = map($_ - 1, @valcolumns_file1);
my @valcolumns_file2 = split(",", $valcolumnstr);
@valcolumns_file2 = map($_ + $offset2, @valcolumns_file2);
print "# columns: ", join("  ", @valcolumns_file2), "  of $file2\n";
@valcolumns_file2 = map($_ - 1, @valcolumns_file2);

my ($triple_val1, $a_ppv1) = tripleidvals($file1, $p1column_file1, \@valcolumns_file1);
my ($triple_val2, $a_ppv2) = tripleidvals($file2, $p1column_file2, \@valcolumns_file2);

print_results( compare_2hashes($triple_val1, $triple_val2 ) );

print "\n";

print_results( compare_2hashes($a_ppv1, $a_ppv2) );


sub print_results{
  my $size1 = shift;
  my $size2 = shift;
  my $n_agree = shift;
  my $n_disagree = shift;
  my $nboth = $n_agree + $n_disagree;
  print "sizes n_1, n_2:            $size1 $size2\n";
  print "n_1only, n_both, n_2only:  ", $size1 - $nboth, "  $nboth  ", $size2 - $nboth, "\n";
  print "n_agree, n_disagree:       $n_agree  $n_disagree\n";
}

  sub tripleidvals{
    my $file = shift;
    my $p1column = shift; 
    my $valcolumns = shift;	# array ref of cols
    open my $fh, "<", "$file";
    my %tripleid_values = ();
    my %aid_p1p2vs = ();
    while (my $line = <$fh>) {
      my @cols = split(" ", $line);
      my $aid = $cols[0];
      next if($aid eq '#');
      my $p1id = $cols[$p1column];
      my $p2id = $cols[$p1column+1];
  
      if (defined $p1id  and  defined $p2id) {
	my $parentids; # = ($p2id gt $p1id)? "$p1id $p2id" : "$p2id $p1id"; #
        my $swap;
	if ($p2id gt $p1id) {
	  $parentids = "$p1id $p2id";
	  $swap = 0;
	} else {
	  $parentids = "$p2id $p1id";
	  $swap = 1;
	}
	my $tripleid = "$aid $parentids";
	if (!exists $tripleid_values{$tripleid}) {
	  $tripleid_values{$tripleid} = '';
	  $aid_p1p2vs{$aid} = "$parentids ";
	}
	my @values = ();
	for my $vcol (@$valcolumns) {
	  my $value = ($vcol < scalar @cols)? $cols[$vcol] : undef;
	  $value = $value // '-'; #
	  push @values, $value;
	}
	if ($swap and  (scalar @values == 5)) {
	  my $t = $values[0]; $values[0] = $values[2]; $values[2] = $t;
	  $t = $values[1]; $values[1] = $values[3]; $values[3] = $t;
	}
	my $vstring = join(" ", @values);
	$tripleid_values{$tripleid} = $vstring;
	$aid_p1p2vs{$aid} .= $vstring;
	# print "$tripleid   $value \n";
      }
    }
    # print scalar keys %tripleid_values, "  ", scalar keys %aid_p1p2vs, "\n";
    close $fh;
    return (\%tripleid_values, \%aid_p1p2vs);
  }

  sub compare_2hashes{ 
    my $h1 = shift;
    my $h2 = shift;
    my ($agree_count, $disagree_count) = (0, 0);
    while (my ($tid, $v) = each %$h1) {
      if (exists $h2->{$tid}) {
	if ($h2->{$tid} eq $v) {
	  $agree_count++;
	} else {
	  $disagree_count++;
	  # print "$tid:  \n", $h2->{$tid}, "\n",   "$v\n"; getc();
	}
      }
    }
    return (scalar keys %$h1, scalar keys %$h2, $agree_count, $disagree_count);
  }
    
