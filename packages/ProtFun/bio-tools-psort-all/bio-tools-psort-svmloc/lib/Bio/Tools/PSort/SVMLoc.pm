package Bio::Tools::PSort::SVMLoc;

use 5.006;
use strict;
use Carp;

require DynaLoader;
require Exporter;
use AutoLoader;

# SVM types
my %SVM_TYPES = ('C-SVC'       => 0,
		 'nu-SVC'      => 1,
		 'one-class'   => 2,
		 'epsilon-SVR' => 3,
		 'nu-SVR'      => 4);
my %SVM_TYPESR = (0 => 'C-SVC',
		  1 => 'nu-SVC',
		  2 => 'one-class',
		  3 => 'epsilon-SVR',
		  4 => 'nu-SVR');

my %KERNEL_TYPES = ('linear'     => 0,
		    'polynomial' => 1,
		    'radial'     => 2,
		    'sigmoid'    => 3);
my %KERNEL_TYPESR = (0 => 'linear',
		     1 => 'polynomial',
		     2 => 'radial',
		     3 => 'sigmoid');

use vars qw(@ISA %EXPORT_TAGS @EXPORT_OK @EXPORT $VERSION);

@ISA = qw(Exporter DynaLoader);

%EXPORT_TAGS = ( 'all' => [ qw( ) ] );
@EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
@EXPORT = qw( );

$VERSION = '0.01';

sub AUTOLOAD {
  my $constname;
  use vars qw($AUTOLOAD);
  ($constname = $AUTOLOAD) =~ s/.*:://;
  croak "& not defined" if $constname eq 'constant';
  my $val = constant($constname, @_ ? $_[0] : 0);
  if ($! != 0) {
    if ($! =~ /Invalid/ || $!{EINVAL}) {
      $AutoLoader::AUTOLOAD = $AUTOLOAD;
      goto &AutoLoader::AUTOLOAD;
    }
    else {
      croak "Your vendor has not defined Bio::Tools::PSort::SVMLoc macro $constname";
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

bootstrap Bio::Tools::PSort::SVMLoc $VERSION;

=head1 NAME

Bio::Tools::PSort::SVMLoc - Perl bindings for the libsvm Support Vector Machine library in Psortb.

=head1 SYNOPSIS

  use Bio::Tools::PSort::SVMLoc;

  # Load the model stored in the file 'sample.model' and Frequent Patterns in 'FreqPatt'
  $svm = new Bio::Tools::PSort::SVMLoc(Model => 'sample.model', FreqPattern => 'FreqPatt');

  # Classify a sequence in $seq, do not include FASTA headers!
  $svm->classify($seq);

=head1 DESCRIPTION

Bio::Tools::PSort::SVMLoc implements a Support Vector Machine for Perl.  Support Vector
Machines provide a method for creating classifcation functions from a set of
labeled training data, from which predictions can be made for subsequent data
sets.

Bio::Tools::PSort::SVMLoc is a stripped down version of Algorithm::SVM design specifically
for PSortb.  If you need functionality such as training of SVM, please use that module
instead.

=head1 CONSTRUCTOR

  # Load an existing SVM.
  $svm = new Bio::Tools::PSort::SVMLoc(Model  => 'sample.model', FreqPattern => 'FreqPatt');

A Bio::Tools::PSort::SVMLoc object must be created loaded from an existing
model file.

An existing SVM is loaded from a file using the Model named parameter and the
frequent patterns set.
The model file should be of the format produced by the svm->train program
(distributed with the libsvm library) or from the $svm->save() method in
Algorithm::SVM.

=head1 METHODS

  $result = $svm->classify($seq);

The classify method is used to classify a set of data according to the
loaded model.  The method accepts a single parameter, which should be
a protein sequence object.  Returns a floating point number
corresponding to the predicted value.

=head1 MAINTAINER

Matthew Laird <matt@brinkman.mbb.sfu.ca>

=head1 SEE ALSO

The libsvm homepage:
http://www.csie.ntu.edu.tw/~cjlin/libsvm/

=head1 ACKNOWLEDGEMENTS

Thanks go out to Fiona Brinkman and the other members of the Simon Fraser
University Brinkman Laboratory for providing me the opportunity to develop
this module.  Additional thanks go to Chih-Jen Lin, one of the libsvm authors,
for being particularly helpful during the development process.

=cut

sub new {
  my ($class, %args) = @_;
  my $self = bless({ }, $class);

  # Ensure we have a valid SVM type.
  $args{Type} = 'C-SVC' if(! exists($args{Type}));
  my $svmtype = $SVM_TYPES{$args{Type}};
  croak("Invalid SVM type: $args{Type}") if(! defined($svmtype));

  # Ensure we have a valid kernel type.
  $args{Kernel} = 'radial' if(! exists($args{Kernel}));
  my $kernel = $KERNEL_TYPES{$args{Kernel}};
  croak("Invalid SVM kernel type: $args{Kernel}") if(! defined($svmtype));

  # Set some defaults.
  my $degree  = exists($args{Degree}) ? $args{Degree} + 0 : 3;
  my $gamma   = exists($args{Gamma}) ? $args{Gamma} + 0 : 0;
  my $coef0   = exists($args{Coef0}) ? $args{Coef0} + 0 : 0;
  my $c       = exists($args{C}) ? $args{C} + 0 : 1;
  my $nu      = exists($args{Nu}) ? $args{Nu} + 0 : 0.5;
  my $epsilon = exists($args{Epsilon}) ? $args{Epsilon} + 0 : 0.1;

  $self->{svm} = _new_svm($svmtype, $kernel, $degree, $gamma, $coef0,
			  $c, $nu, $epsilon);

  # Load the model if one was specified.
  my $model = $args{Model};
  croak("Model file not found or bad permissions: $model")
      if((! -r $model) || (! -f $model));

  # Load the model.
  croak("Error loading the model: $model")
      unless($self->loadModel($model));

  my $freqpatterns = $args{FreqPatt};

  # Load the frequent patterns if we have any
  if($freqpatterns) {
      croak("Error loading the frequent patterns: $freqpatterns")
	  unless($self->loadFreqPattern($freqpatterns));
  }

  # Ensure that the model loaded correctly.
  croak("Error loading model file: $model") if(! $self->{svm});

  return $self;
}

sub loadModel {
  my ($self, $file) = @_;

  croak("Can't load model because no filename provided") if(! $file);

  return _loadModel($self->{svm}, $file);
}

sub loadFreqPattern {
  my ($self, $file) = @_;

  croak("Can't load model because no filename provided") if(! $file);

  croak("Freq Patterns file not found or bad permissions: $file: $!")
      if((! -r $file) || (! -f $file));

  return _loadFreqPattern($self->{svm}, $file);
}

sub classify {
  my ($self, $x) = @_;

  return _classify($self->{svm}, $x);
 }

sub predict {
  my ($self, $x) = @_;

  # Check if we got a dataset object.
  croak("Not an  Bio::Tools::PSort::SVMLoc::DataSet") if(ref($x) ne "Bio::Tools::PSort::SVMLoc::DataSet");

  return _predict($self->{svm}, $x);
 }

1;

__END__
