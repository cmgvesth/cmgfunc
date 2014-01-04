package Algorithm::HMM::Report;

use strict;
use 5.006;
use Carp;

use Algorithm::HMM::Hit;

=head1 NAME

Algorithm::HMM::Report - Report object for the Algorithm::HMM package.

=head1 DESCRIPTION

The Algorithm::HMM::Report object is returned by the Algorithm::HMM search method.  It contains the sets of global and domain hits returned by the model.

=head1 SYNOPSIS

  # Search the HMM model.  An Algorithm::HMM::Report object is returned.
  my $rep = $hmm->search("AAIELKBPOWELKQJPASDLKJIGE");

  # Get all the global hits from the report object.
  my @ghits = $rep->global_hits();

  # Get all the domain hits from the report object.
  my @dhits = $rep->domain_hits();

=head1 METHODS

  my @ghits = $rep->global_hits();
  my @dhits = $rep->domain_hits();

The global_hits() and domain_hits() methods return lists of Algorithm::HMM::Hit::Global and Algorithm::HMM::Hit::Domain objects respectively.  A global hit is a match across the entire model, while a domain hit is a match on a smaller subsequence of the model.

=head1 AUTHOR

The Algorithm::HMM package was originally written by Cory Spencer <cspencer@sfu.ca> of the Simon Fraser University Brinkman Laboratory.  It is currently maintained by Matthew Laird <matt@brinkman.mbb.sfu.ca>.

=head1 SEE ALSO

Algorithm::HMM, Algorithm::HMM::Hit::Global, Algorithm::HMM::Hit::Domain

=cut

sub new {
  my ($class, %args) = @_;
  my $self = bless({ }, $class);

  # Get the list of global hits and make sure it's an array ref.
  croak("Global hits must be array ref")
    if(defined($args{Global}) && (ref($args{Global}) ne "ARRAY"));
  $self->{global} = $args{Global} || [ ];

  # Get the list of global hits and make sure it's an array ref.
  croak("Domain hits must be array ref")
    if(defined($args{Domain}) && (ref($args{Domain}) ne "ARRAY"));
  $self->{domain} = $args{Domain} || [ ];

  return $self;
}

sub global_hits {
  my $self = shift;
  my @hits;

  push(@hits, $self->_getGlobalHit($_)) for(0..($self->_numGlobalHits() - 1));
  return @hits;
}

sub domain_hits {
  my $self = shift;
  my @hits;

  push(@hits, $self->_getDomainHit($_)) for(0..($self->_numDomainHits() - 1));
  return @hits;
}

1;

__END__
