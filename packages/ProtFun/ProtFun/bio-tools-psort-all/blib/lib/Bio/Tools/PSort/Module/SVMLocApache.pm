package Bio::Tools::PSort::Module::SVMLocApache;

use Bio::Tools::PSort::Constants qw(:all);
use Bio::Tools::PSort::Module::AnalysisI;
use Bio::Tools::PSort::Report::Result;
use LWP;

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
  my ($url, $module) = $self->_rearrange([qw(URL MODULE)], @args);

  $self->{url} = $url ||
      $self->throw("URL seems to be blank, where are we pointing?");

  $self->{module} = $module ||
      $self->throw("Module seems to be blank, what are we fetching?");

  $self->{browser} = LWP::UserAgent->new ||
      $self->throw("We weren't able to make an LWP Agent, this is bad.");


  return $self;
}

sub run {
  my ($self, $seq) = @_;


  $self->throw("Not a Bio::Seq object")
    if((! ref($seq)) && (! $seq->isa("Bio::Seq")));

  $self-throw("No sequence in Bio::Seq object")
      unless($seq->seq);

  my $response = $self->{browser}->post ( $self->{url},
					  [ 'Module' => $self->{module},
					    'Seq' => $seq->seq 
					  ] 
					);

  $self->throw("Error, wrong content type from SVM Module $self->{module} at $self->{url}, received " . $response->content_type)
      unless $response->content_type eq 'text/plain';

#  my $loc = $LOCS{'Unknown'};
  my $loc = -1;
  # Make sure we got a number, if this files $1 will be ''
  $response->content =~ /$self->{module}:\s+(\d+)/m;
  my $loc = $1;
  # Ensure the result is a number and that the previous check did cause it to be ''
  $loc = ($loc =~ /\d+/ ? $loc : 10);
#  $loc = $1 || 10;

  return new Bio::Tools::PSort::Report::Result(-details => [ ],
					       -score   => 0,
					       -loc     => $LOCSR{$loc});
}

1;
