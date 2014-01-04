package Bio::Tools::PSort::Module::SVMLoc;

use Bio::Tools::PSort::Constants qw(:all);
use Bio::Tools::PSort::Module::AnalysisI;
use Bio::Tools::PSort::Report::Result;
use Bio::Tools::PSort::SVMLoc;

use vars qw(@ISA);
use strict;

@ISA = qw(Bio::Tools::PSort::Module::AnalysisI);


my %LOCMAP = (Cytoplasmic  => CYTOPLASM,
	      CytoplasmicMembrane => CYTOPLASMICMEMBRANE,
	      Periplasmic   => PERIPLASMIC,
	      OuterMembrane => OUTERMEMBRANE, 
	      Extracellular => EXTRACELLULAR,
	      Mitochondrial => MITOCHONDRIAL,
	      Nuclear       => NUCLEAR,
	       Unknown => UNKNOWN,
	      Membrane => MEMBRANE,
	      Cellwall => CELLWALL, 
	      Other => OTHER	   
	 	    
);

my %LOCS  = (Cytoplasmic   => 0, CytoplasmicMembrane => 1, Periplasmic   => 2,
	     OuterMembrane => 3, Extracellular => 4, Mitochondrial => 5,
	     Nuclear       => 6, Other       => 7,
	     Membrane => 8, Cellwall => 9, Unknown => 10 );

my %LOCSR = (0 => 'Cytoplasmic',   1 => 'CytoplasmicMembrane', 2 => 'Periplasmic',
	     3 => 'OuterMembrane', 4 => 'Extracellular', 5 => 'Mitochondrial',
	     6 => 'Nuclear',       7 => 'Other', 
	     8 => 'Membrane', 9 => 'Cellwall', 10 => 'Unknown');

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

##receive the SVM_MODEL (based on training data) and the Frequent patterns file
  my ($model, $patterns) = $self->_rearrange([qw(MODEL PATTERNS)], @args);

  $self->{subloc} = new Bio::Tools::PSort::SVMLoc(Model => $model, FreqPatt => $patterns) ||
    $self->throw("Couldn't create new Bio::Tools::PSort::SVMLoc object");

  return $self;
}

sub run {
  my ($self, $seq) = @_;


  $self->throw("Not a Bio::Seq object")
    if((! ref($seq)) && (! $seq->isa("Bio::Seq")));

  $self-throw("No sequence in Bio::Seq object")
      unless($seq->seq);

  my $loc = $self->{subloc}->classify($seq->seq);
 
  return new Bio::Tools::PSort::Report::Result(-details => [ ],
					       -score   => 0,
					       -loc     => $LOCSR{$loc});
}

1;
