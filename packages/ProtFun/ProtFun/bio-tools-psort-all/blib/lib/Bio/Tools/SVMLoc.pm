package Bio::Tools::SVMLoc;

use Algorithm::SVM::DataSet;
use Algorithm::SVM;

use Bio::Root::Root;
use Carp;

# SVMLoc localizations.
my %LOCS  = (Cytoplasmic   => 0, CytoplasmicMembrane => 1, Periplasmic   => 2,
	     OuterMembrane => 3, Extracellular => 4, Mitochondrial => 5,
	     Nuclear       => 6, Other       => 7,
	     Membrane => 8, Cellwall => 9, Unknown => 10 );

my %LOCSR = (0 => 'Cytoplasmic',   1 => 'CytoplasmicMembrane', 2 => 'Periplasmic',
	     3 => 'OuterMembrane', 4 => 'Extracellular', 5 => 'Mitochondrial',
	     6 => 'Nuclear',       7 => 'Other', 
	     8 => 'Membrane', 9 => 'Cellwall', 10 => 'Unknown');

use strict;

our @ISA = qw(Bio::Root::Root);
our $VERSION = '0.03';

# List of all the amino acid symbols.
my @AA = qw(V L I M F W Y G A P S T C H R K Q E N D);

sub new {

  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  # Get the model filename.
  my ($model) = $self->_rearrange([qw(MODEL)], @args);

  # Load model from a file if specified.
  if($model) {
      $self->{svm} = new Algorithm::SVM(Model => $model);
    $self->throw("SVMLoc SVM initialization failed") if(! $self->{svm});
  }

  return $self;
}

sub classify {


  my ($self,  $seq, $fre_patterns) = @_;
  my ($ds, $loc);

  # Ensure the frequent patterns is loaded
  $self->throw("No frequent patterns loaded") if(!(ref($fre_patterns) eq 'ARRAY'));
  my @patterns = @$fre_patterns;
 
  # Ensure we have a model loaded.
  $self->throw("No model is currently loaded") if(! $self->{svm});

  # Make sure that we received a sequence object.
  $self->throw("Not a Bio::Seq object")
      if((! ref($seq)) || (! $seq->isa("Bio::SeqI")));

  # Create a new dataset to use in our test against the SVM model
  $ds = new Algorithm::SVM::DataSet(Label => 1);

  #value is set to 1 if the patterns are present
  my $value = 1;
 
  #iterate through all of the sequences, and look for whether each of the patterns is present in the sequence
  my $index;
  if($seq->seq) {
      for $index (0 .. $#patterns) {
	  if((index $seq->seq, $patterns[$index]) != -1) {
	      $ds->attribute($index+1, $value);
	  }
      }
  }
      
  # Make the SVM prediction.
  $loc = $self->{svm}->predict($ds);

  return $LOCSR{$loc};

}

sub train {

  # $data is a reference to an array containing the vectors from the train.dat file
  # $SCL is a string, the localization
  # $kernel is what SVM kernel to use
  my ($self, $data, $SCL, $opts) = @_;

  my @training = @$data;
  my $ds;
  my @tset = ();

  # Check if we're overriding the kernel
  my $kernel = $opts->{Kernel} || 'linear';

  $self->{svm} = new Algorithm::SVM(Type => 'C-SVC', Kernel => $kernel, Gamma => 1);
  $self->throw("Error initializing SVM") if(! $self->{svm});

  foreach my $line (@training) {
      chomp $line;

      if($line =~ /^1/) {
	  $ds = new Algorithm::SVM::DataSet(Label => $LOCS{$SCL});
      } elsif($line =~ /^-1/) {
	  $ds = new Algorithm::SVM::DataSet(Label => $LOCS{'Unknown'});
      }

      my @vector = split ' ', $line;
      for my $i (1 .. $#vector) {
	  my ($index) = split ':', $vector[$i];
	  $ds->attribute($index, 1);
      }
      
      push(@tset, $ds);
	  
  }

  if($opts->{Optimize}) {
      $self->{svm}->train(@tset);

      my ($high_g, $high_c, $high_score);

      $high_score = 0;

      my $start_g = $opts->{GammaStart} || 1;
      my $end_g   = $opts->{GammaEnd}   || 100;
      my $inc_g   = $opts->{GammaInc}   || 1;

      my $start_c = $opts->{CStart} || 1;
      my $end_c   = $opts->{CEnd}   || 50;
      my $inc_c   = $opts->{CInc}   || 1;

      for(my $i = $start_g; $i <= $end_g; $i += $inc_g) {
	  for(my $j = $start_c; $j <= $end_c; $j += $inc_c) {
	      $self->{svm}->gamma($i);
	      $self->{svm}->C($j);

	      $self->{svm}->retrain();

	      my $score = $self->{svm}->validate(5);
	      if($score > $high_score) {
		  $self->debug("Accuracy: $score% (Gamma = $i, C = $j)\n");

		  $high_score = $score;
		  $high_g = $i;
		  $high_c = $j;
	      }
	  }
      }

      $self->{svm}->gamma($high_g);
      $self->{svm}->C($high_c);
      $self->{svm}->retrain();
  } else {
      return $self->{svm}->train(@tset);
  }

}


sub validate {

    my ($self) = @_;
    print "performing 5 fold validation....\n";
    my $accuracy = $self->{svm}->validate(5);
    print "ACCURACY =  $accuracy\n";
}


sub save {
  my ($self, $filename) = @_;

  return $self->{svm}->save($filename);
}

sub load {
  my ($self, $filename) = @_;

  return $self->{svm}->load($filename);
}


=head1 NAME

Bio::Tools::SVMLoc - Perl implementation of the SVM based protein subcellular
localization method.

=head1 SYNOPSIS

  use Bio::Tools::SVMLoc;
  use Bio::SeqIO;

  # Load a previously trained model from a file.
  $svmloc = new Bio::Tools::SVMLoc(-model  => 'svmloc.model');

  # Classify a Bio::Seq object.
  $loc = $svmloc->classify($seq);

  # Train SVMLoc on a new dataset.
  $accuracy = $svmloc->train(\@vectors, $localization, 
                             { Kernel  => 'linear',
                               Optimize => 1 });

  # Save the model to a file.
  $svmloc->save('svmloc.model.new');

  # Load a model from a file.
  $svmloc->load('svmloc.model.new');

=head1 DESCRIPTION

Bio::Tools::SVMloc is an implementation of the SVMLoc protein subcellular
localization method, which predicts localizations based on amino acid
composition.  The method was originally outlined in "Support Vector
Machine Approach for Protein Subcellular Localization Prediction" paper by
Sujun Hua and Zhirong Sun.

=head1 CONSTRUCTOR

   $svmloc = new Bio::Tools::SVMLoc(-model  => 'svmloc.model');

The SVMLoc constructor accepts the name of an existing model file.

=head1 METHODS

   $loc = $sl->classify($seq, \@frequent_patterns);

The classify method accepts a Bio::Seq object and an array reference of
the known frequent patterns as arguments and returns a string containing
the predicted localization of the protein.  The return value will be one
of: Cytoplasmic, CytoplasmicMembrane, Periplasmic, OuterMembrane, Secreted,
Mitochondrial, Nuclear or Unknown.

  $accuracy = $svmloc->train(\@vectors, $localization, 
                             { Kernel  => 'linear',
                               Optimize => 1 });

The train method allows the creation of a new SVM model based on a
set of user defined sequences.  The first parameter an array reference
containing the list of vectors to train the SVM on, those in the 
positive set will be prefixed by 1 while those in the negative set
will be prefixed by -1.  The second parameter is a string containing
the localization that is being targeted with this SVM.

The resulting model is unlikely to provide optimal results without some
degree of fine tuning to the underlying Support Vector Machine.  The
third (optional) parameter to the train method is a hashref which allows
one to pass options directly to the SVM.  The set of valid keys are as
follows:

  Kernel     - The Support Vector Machine kernel to use, possible values
               are linear, polynomial, radial, and sigmoid.
  
  Optimize   - If true, the module will attempt to search for the best
               values of the gamma and C parameters passed to the SVM.
               For any new model, this should be run at least once.

               By default, the SVMLoc module will try all values for gamma
               between 1 and 100 (increments of 1) and C between 1 and 50
              (increments of one) until the optimal values of each are
               located.

               Note that this **WILL** take quite some time, but only
               needs to be done once per model.  The model should be
               saved with the save method after training so that all
               optimal parameter values will be retained.

  GammaStart - Allows specification of the start value when searching for
               the optimal gamma value.  (Default 1)

  GammaEnd   - Allows specification of the end value when searching for
               the optimal gamma value.  (Default 100)

  GammaInc   - Allows specification of the size of the increment when
               searching for the optimal gamma value.  (Default 1)

  CStart      - Allows specification of the start value when searching for
               the optimal C value.  (Default 1)

  CEnd        - Allows specification of the end value when searching for
                the optimal C value.  (Default 50)

  CInc        - Allows specification of the size of the increment when
                searching for the optimal C value.  (Default 1)

  $svmloc->save('svmloc.model.new');

Saves the currently loaded model to the specified file.  Returns true on
success, false on failure.

  $svmloc->load('svmloc.model.new');

Loads an existing model from file.  Returns true on success, false on
failure.


=head1 AUTHOR

Fiona Brinkman, Cory Spencer, Brinkman Lab, Simon Fraser University <psort-mail@sfu.ca>

=head1 MAINTAINER

Matthew Laird <lairdm@sfu.ca>

=head1 SEE ALSO

Algorithm::SVM and the original SVMLoc paper titled "Support Vector Machine
Approach for Protein Subcellular Localization Prediction" by
Sujun Hua and Zhirong Sun

=head1 ACKNOWLEGEMENTS

Chih-Jen Lin, one of the authors of the libsvm package upon which
Algorithm::SVM was based, was invaluable in creating this version of SVMLoc
Thanks also to Sujun Hua and Zhirong Sun, the authors of the original SVMLoc.
Last but not least, thanks to Fiona Brinkman, Jennifer Gardy and the other
members of the Simon Fraser University Brinkman laboratory.

=cut

1;
