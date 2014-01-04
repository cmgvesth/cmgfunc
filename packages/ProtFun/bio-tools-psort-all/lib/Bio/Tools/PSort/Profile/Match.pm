package Bio::Tools::PSort::Profile::Match;

use Bio::Root::Root;
use Carp;

use strict;

our @ISA = qw(Bio::Root::Root);
our $VERSION = '0.01';

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  # Get the parameters we've been passed.
  my ($id,$loc,$comment,$start,$end) =
    $self->_rearrange([qw(PROFILEID LOCALIZATION COMMENT START END)], @args);

  $self->{localization} = $loc;
  $self->{profile_id}     = $id;
  $self->{comment}      = $comment;
  $self->{match_start}  = $start;
  $self->{match_end}    = $end;

  return $self;
}

sub profile_id {
  my $self = shift;

  return (@_) ? $self->{profile_id} = shift : $self->{profile_id};
}

sub localization {
  my $self = shift;

  return (@_) ? $self->{localization} = shift : $self->{localization};
}

sub comment {
  my $self = shift;

  return (@_) ? $self->{comment} = shift : $self->{comment};
}

sub start {
  my $self = shift;

  return (@_) ? $self->{match_start} = shift : $self->{match_start};
}

sub end {
  my $self = shift;

  return (@_) ? $self->{match_end} = shift : $self->{match_end};
}

=head1 NAME

  Bio::Tools::Motif::Match - Represents a match as reported by the
  Bio::Tools::Motif module.

=head1 SYNOPSIS

  use Bio::Tools::Motif::Match;

  # Create a new object.
  $match = new Bio::Tools::Motif::Match(-motifid => 'P07110',
                                        -localization => 'OuterMembrane',
                                        -comment => 'FIMBRIAL_USHER Pattern',
                                        -start => 314,
                                        -end => 325);

  # Accessor methods.
  $val = $match->motif_id();
  $val = $match->localization();
  $val = $match->comment();
  $val = $match->start();
  $val = $match->end();

=head1 AUTHOR

Cory Spencer <cspencer@sfu.ca>

=head1 SEE ALSO

  Bio::Tools::Motif

=head1 ACKNOWLEGEMENTS

Thanks go out to Fiona Brinkman, Jennifer Gardy and the other members of the
Simon Fraser University Brinkman laboratory.

=cut

1;
