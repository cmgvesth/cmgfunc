package Algorithm::HMM::Hit::Domain;

use strict;
use 5.006;
use Carp;

use Algorithm::HMM::Hit;

our @ISA = qw(Algorithm::HMM::Hit);

=head1 NAME

Algorithm::HMM::Hit::Domain - Object encapsulating information about a Hidden Markov Model domain hit.

=head1 DESCRIPTION

The Algorithm::HMM::Hit::Domain object encapsulates information about a domain hit in the Hidden Markov Model.  Unlike global hits which attempt to match the provided sequence against the entire model, a domain hit indicates a match against a smaller subsequence in the model.

=head1 METHODS

  $dhit->pvalue()
  $dhit->evalue()
  $dhit->score()

Returns the pvalue, evalue and score values respectively for the hit.

  $dhit->seq_from()
  $dhit->seq_to()
  $dhit->seq_length()

Returns the starting and ends points in the sequence that matched the model.  seq_length returns the length of the hit.

  $dhit->hmm_from()
  $dhit->hmm_to()

Returns the starting and end points in the model.

  $dhit->parent_pvalue()
  $dhit->parent_evalue()
  $dhit->parent_score()

Returns the parent pvalue, evalue and score values respectively.

=head1 AUTHOR

The Algorithm::HMM package was originally written by Cory Spencer <cspencer@sfu.ca> of the Simon Fraser University Brinkman Laboratory.  It is currently maintained by Matthew Laird <matt@brinkman.mbb.sfu.ca>.

=head1 SEE ALSO

Algorithm::HMM, Algorithm::HMM::Hit::Global, Algorithm::HMM::Report

=cut

sub new {
  my ($class, %args) = @_;
  my $self = $class->SUPER::new(%args);


  return $self;
}

sub seq_from {
  my $self = shift;

  return (@_) ? $self->_setSeqFrom(shift(@_)) : $self->_getSeqFrom();
}

sub seq_to {
  my $self = shift;

  return (@_) ? $self->_setSeqTo(shift(@_)) : $self->_getSeqTo();
}

sub seq_length {
  my $self = shift;

  return (@_) ? $self->_setSeqLength(shift(@_)) : $self->_getSeqLength();
}

sub hmm_from {
  my $self = shift;

  return (@_) ? $self->_setHMMFrom(shift(@_)) : $self->_getHMMFrom();
}

sub hmm_to {
  my $self = shift;

  return (@_) ? $self->_setHMMTo(shift(@_)) : $self->_getHMMTo();
}

sub parent_pvalue {
  my $self = shift;

  return (@_) ? $self->_setParentPValue(shift(@_)) : $self->_getParentPValue();
}

sub parent_evalue {
  my $self = shift;

  return (@_) ? $self->_setParentEValue(shift(@_)) : $self->_getParentEValue();
}

sub parent_score {
  my $self = shift;

  return (@_) ? $self->_setParentScore(shift(@_)) : $self->_getParentScore();
}

1;

__END__
