package Algorithm::HMM;

use strict;
use 5.006;
use Carp;

require DynaLoader;
require Exporter;
use AutoLoader;

our @ISA = qw(Exporter DynaLoader);

our %EXPORT_TAGS = ( 'all' => [ qw( ) ] );
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
our @EXPORT = qw( );

our $VERSION = '0.02';

sub AUTOLOAD {
  my $constname;
  our $AUTOLOAD;
  ($constname = $AUTOLOAD) =~ s/.*:://;
  croak "& not defined" if $constname eq 'constant';
  my $val = constant($constname, @_ ? $_[0] : 0);
  if ($! != 0) {
    if ($! =~ /Invalid/ || $!{EINVAL}) {
      $AutoLoader::AUTOLOAD = $AUTOLOAD;
      goto &AutoLoader::AUTOLOAD;
    }
    else {
      croak "Your vendor has not defined Algorithm::HMM macro $constname";
    }
  }
  {
    no strict 'refs';
    # Fixed between 5.005_53 and 5.005_61
    if ($] >= 5.00561) {
      *$AUTOLOAD = sub () { $val };
    }
    else {
      *$AUTOLOAD = sub { $val };
    }
  }
  goto &$AUTOLOAD;
}

bootstrap Algorithm::HMM $VERSION;

use Algorithm::HMM::Hit::Domain;
use Algorithm::HMM::Hit::Global;
use Algorithm::HMM::Report;

=head1 NAME

Algorithm::HMM - Perl bindings for a Hidden Markov Model library.

=head1 SYNOPSIS

  use Algorithm::HMM;

  # Load a Hidden Markov Model.
  my $hmm = new Algorithm::HMM(Model => 'sample.hmm');

  # Run a sequence through the loaded model.
  my $rep = $hmm->search("AAIELKBPOWELKQJPASDLKJIGE");

  # Get all the global hits from the search.
  my @ghits = $rep->global_hits();

  # Get all the domain hits from the search.
  my @dhits = $rep->domain_hits();

  # Display information about the domain hits.
  for (@dhits) {
    print("pvalue = " . $_->pvalue . "\n");
    print("evalue = " . $_->evalue . "\n");
    print("score  = " . $_->score  . "\n");
  }

  # Save the model to a different file.
  $hmm->save("sample.hmm.0");

  # Load the saved model again.
  $hmm->load("sample.hmm.1")

=cut

sub new {
  my ($class, %args) = @_;
  my $self = bless({ }, $class);

  my $model = $args{Model} || '';
  my $df    = ($args{DoForward} || 0) + 0;
  my $n2    = ($args{DoNull2} || 0) + 0;

  $self->{hmm} = hmm_new($model, $df, $n2);
  croak("Couldn't create new HMM") if(! $self->{hmm});

  return $self;
}

sub load {
  my ($self, $file) = @_;

  croak("Can't load model because no filename provided") if(! $file);

  return $self->{hmm}->_load($file);
}

sub save {
  my ($self, $file, %args) = @_;
  my ($append, $binary);

  croak("Can't save model because no filename provided") if(! $file);
  croak("No model currently loaded") if(! $self->{hmm}->_modelLoaded());

  $append = $args{Append} + 0 || 0;
  $binary = $args{Binary} + 0 || 0;

  return $self->{hmm}->_save($file, $append, $binary);
}

sub search {
  my ($self, $seq) = @_;
  my ($rep, @global, @domain);

  $rep = $self->{hmm}->_search($seq);
  croak("Failed to create Algorithm::HMM::Report") if(! defined($rep));

  return $rep;
}

sub forward {
  my $self = shift;

  (@_) ? _setDoForward($self->{hmm},shift(@_)+0) : _getDoForward($self->{hmm});
}

sub null2 {
  my $self = shift;

  (@_) ? _setDoNull2($self->{hmm},shift(@_)+0) : _getDoNull2($self->{hmm});
}

sub train {
  my ($self, @seqs) = @_;
  my $max;

  # Determine the maxium sequence length.
  for(@seqs) { $max = length($_) if($max < length($_)) }

  # Pad the sequences to the same length.
  $_ = $_ . "-" x ($max - length($_)) for(@seqs);

  return _train($self->{hmm}, @seqs);
}

1;

__END__
