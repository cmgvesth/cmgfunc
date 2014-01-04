package Bio::Tools::PSort::SVMLoc::DataSet;

use 5.006;
use strict;
use Carp;

use Bio::Tools::PSort::SVMLoc;

=head1 NAME

Bio::Tools::PSort::SVMLoc::DataSet - A DataSet object for the Bio::Tools::PSort::SVMLoc Support
Vector Machine.

=head1 SYNOPSIS

  use Bio::Tools::PSort::SVMLoc::DataSet;

  # Create a new dataset.
  $ds = new Bio::Tools::PSort::SVMLoc::DataSet(Label => 1,
                                    Data  => [ 0.12, 0.25, 0.33, 0.98 ]);

  # Retrieve/set the label.
  $label = $ds->label();
  $ds->label(1976);

  # Retrieve/set the attribute with an index of 0.
  $attr = $ds->attribute(0);
  $ds->attribute(0, 0.2621);

=head1 DESCRIPTION

Bio::Tools::PSort::SVMLoc::DataSet is a representation of the datasets passed to
Bio::Tools::PSort::SVMLoc object for training or classification.  Each dataset has
an associated label, which classifies it as being part of a specific group.
A dataset object also has one or more key/value pairs corresponding to
the attributes that will be used for classification.

=head1 CONSTRUCTORS

 $ds = new Bio::Tools::PSort::SVMLoc::DataSet(Label => 1,
                                   Data  => [ 0.12, 0.25, 0.33, 0.98 ]);

The Bio::Tools::PSort::SVMLoc::DataSet constructor accepts two optional named 
parameters: Label and Data.  Label is used to set the class to which the
dataset belongs, and Data is used to set any initial values.  Data
should be an arrayref of numerical values.  Each value in the arrayref
is assumed to have a key corresponding to its index in the array.

  ie) In the above example, 0.12 has a key of 0, 0.25 has a key of 1,
      0.33 has a key of 2, etc.

=head1 METHODS

  $label = $ds->label();
  $ds->label(1976);

The label method is used to set or retrieve the DataSets label value.
Parameters and return values should be numeric values.

  $attr = $ds->attribute(0);
  $ds->attribute(0, 0.2621);


The attribute method is used to set dataset attribute values.  If a single
value is provided, the method will return the corresponding value.  If
two value are provided, the method will set the first parameter to the
value of the second.

=head1 MAINTAINER

Matthew Laird <matt@brinkman.mbb.sfu.ca>

=head1 SEE ALSO

Bio::Tools::PSort::SVMLoc

=cut


sub new {
  my ($class, %args) = @_;

 # Do some quick error checking on the values we've been passed.
  croak("No label specified for DataSet") if(! exists($args{Label}));
  my $self = _new_dataset($args{Label} + 0);

  if(exists($args{Data})) {
    croak("Data must be an array ref") if(ref($args{Data}) ne "ARRAY");
    for(my $i = 0; $i < @{$args{Data}}; $i++) {
      $self->attribute($i, (@{$args{Data}})[$i] + 0);
    }
  }

  return $self;
}

sub label {
  my ($self, $label) = @_;

  return (defined($label)) ? _setLabel($self, $label + 0) : _getLabel($self);
}

sub attribute {
  my ($self, $key, $val) = @_;

  croak("No key specified") if(! defined($key));

  return(defined($val)) ? _setAttribute($self, int($key), $val + 0) :  _getAttribute($self, int($key));
}

1;

__END__
