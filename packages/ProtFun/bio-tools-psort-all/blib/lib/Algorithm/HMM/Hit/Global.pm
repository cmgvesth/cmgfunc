package Algorithm::HMM::Hit::Global;

use strict;
use 5.006;
use Carp;

use Algorithm::HMM::Hit;

our @ISA = qw(Algorithm::HMM::Hit);

=head1 NAME

Algorithm::HMM::Hit::Global - Object encapsulating information about a Hidden Markov Model global hit.

=head1 DESCRIPTION

The Algorithm::HMM::Hit::Global object encapsulates information about a global hit in the Hidden Markov Model.  Unlike domain hits are matches against smaller subsequences in the model, a global hit indicates a match against the entire model.

=head1 METHODS

  $dhit->pvalue()
  $dhit->evalue()
  $dhit->score()

Returns the pvalue, evalue and score values respectively for the hit.

=head1 AUTHOR

The Algorithm::HMM package was originally written by Cory Spencer <cspencer@sfu.ca> of the Simon Fraser University Brinkman Laboratory.  It is currently maintained by Matthew Laird <matt@brinkman.mbb.sfu.ca>.

=head1 SEE ALSO

Algorithm::HMM, Algorithm::HMM::Hit::Domain, Algorithm::HMM::Report

=cut

sub new {
  my ($class, %args) = @_;
  my $self = $class->SUPER::new(%args);


  return $self;
}

sub DESTROY {

}

1;

__END__
