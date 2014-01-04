package Algorithm::HMM::Hit;

use strict;
use 5.006;
use Carp;

use Algorithm::HMM;

=head1 NAME

Algorithm::HMM::Hit - Abstract class implementing methods shared by the Algorithm::HMM::Hit::Domain and Algorithm::HMM::Hit::Global classes.

=head1 AUTHOR

The Algorithm::HMM package was originally written by Cory Spencer <cspencer@sfu.ca> of the Simon Fraser University Brinkman Laboratory.  It is currently maintained by Matthew Laird <matt@brinkman.mbb.sfu.ca>.

=head1 SEE ALSO

Algorithm::HMM, Algorithm::HMM::Hit::Global, Algorithm::HMM::Domain, Algorithm::HMM::Report

=cut

use strict;

our $VERSION = '0.01';

sub new {
  my $class = shift;

  croak("Can't instantiate abstract class: $class");
}

sub pvalue {
  my $self = shift;

  return $self->_getPValue();
}

sub evalue {
  my $self = shift;

  return $self->_getEValue();
}

sub score {
  my $self = shift;

  return $self->_getScore();
}

sub DESTROY {

}

1;

__END__
