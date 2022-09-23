package PedigreesDAG;
use strict;
use warnings;
use Moose;
#use Mouse;
use namespace::autoclean;
use Carp;
use Scalar::Util qw (looks_like_number );
use List::Util qw ( min max sum );
use POSIX qw ( floor ceil );
use Genotypes;
use Node;
#use DAG;

has pedigree_filename => (
			  isa => 'Str',
			  is => 'ro',
			  required => 1,
			 );

# has accid_parents => (
# 		      isa => 'HashRef', # keys are accession ids, values: [matid, patid]
# 		      is => 'rw',
# 		      default => sub{ {} },
# 		     );

has id_node => (
		isa => 'HashRef[Maybe[Str]]', # keys are accession ids, values are Node objects.
		is => 'rw',
		default => sub{ {} },
	       );

has childless_ids => (
		      isa => 'HashRef[Maybe[Str]]', #  key: id (string), value: Node Obj
		      is => 'rw',
		      default => sub{ {} },
		     );

# has parentless_ids => (
# 		       isa => 'HashRef[Maybe[Str]]', #  key: id (string), value: Node Obj
# 		       is => 'rw',
# 		       default => sub{ {} },
# 		      );


sub BUILD{
  my $self = shift;
  my $pedigree_filename = $self->pedigree_filename();

  #  my $accid_parentalids = {};
  open my $fh, "<", "$pedigree_filename" or die "couldn't open $pedigree_filename for reading.\n";
#  my $first_line = <$fh>;
#  die "pedigree file should have 'Accession' at start of first line.\n" if(! ($first_line =~ /^Accession/));
  #  my $the_dag = DAG->new();
  while (my $line = <$fh>) {
    my @cols = split(" ", $line);
    my ($accid, $matid, $patid) = @cols[-3, -2, -1];
    next if($accid eq 'NA');

    my $matnode = ($matid eq 'NA')? undef : ( $self->id_node()->{$matid} // Node->new({id => $matid}) );
    my $patnode = ($patid eq 'NA')? undef : ( $self->id_node()->{$patid} // Node->new({id => $patid}) );
    $matnode->add_offspring( $accid ) if(defined $matnode);
    $patnode->add_offspring( $accid ) if(defined $patnode);

    my $anode = $self->id_node()->{$accid} // Node->new({id => $accid, female_parent => $matnode, male_parent => $patnode});
    $self->add_node($anode);

  }

  while (my ($id, $anode) = each  %{$self->id_node()}) {
    if ($anode->offspring() eq '') {
      $self->childless_ids()->{$id} = 1;
    }
  }

 
  my ($acyclic, $cyclic_id_string) = $self->is_it_acyclic();
  my $info_string = "# Is it acyclic:  $acyclic \n";
  $info_string .= "# Childless accessions with ancestral cyclicities: $cyclic_id_string \n" if($acyclic == 0);
  $info_string .= "# N ids: ", scalar keys %{$self->id_node()}, " N with no offspring: ", scalar keys %{$self->childless_ids()}, "\n";
}

sub add_node{
  my $self = shift;
  my $node = shift;
  my $node_id = $node->id();
  $self->id_node()->{$node_id} = $node;
}

sub parents{			# return ids of parents, if defined.
  my $self = shift;
  my $accid = shift;
  if (defined (my $the_node = $self->id_node()->{$accid})){ # // undef;
    my $fpid = (defined $the_node->female_parent())? $the_node->female_parent()->id() : undef;
    my $mpid = (defined $the_node->male_parent())? $the_node->male_parent()->id() : undef;
    return ($fpid, $mpid);
  } else {
    return (undef, undef);
  }
}

sub as_string{
  my $self = shift;
  my $string = '';
  my @ids = sort keys %{$self->id_node()};
  for my $id (@ids) {
    my $node = $self->id_node()->{$id};
    my ($nwck, $mind, $maxd) = $node->as_newick();
    my @offspring_array = split(" ", $node->offspring());
    $string .= "$id  $mind $maxd   $nwck  " . scalar @offspring_array . "\n";
  }
  return $string;
}

sub is_it_acyclic{
  my $self = shift;
  my $cyclic_nodes = '';	#
  my $acyclicity = 1;
  for my $anode ( map($self->id_node()->{$_}, sort keys %{$self->id_node()} ) ) {
    next if($anode->offspring() ne ''); # skip if has offspring
    my $node_acyclicity = $anode->ancestors_acyclic( {} ); # is directed graph of node and ancestors acyclic
    if ($node_acyclicity == 0) {
      $acyclicity = 0;
      $cyclic_nodes .= $anode->id() . ' ';
    }
  }
  return ($acyclicity, $cyclic_nodes);
}

1;
