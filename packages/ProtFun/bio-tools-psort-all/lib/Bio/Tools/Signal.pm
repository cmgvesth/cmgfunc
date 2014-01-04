package Bio::Tools::Signal;

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# OVERVIEW

# PSORT-B is described in Gardy, J.L. et al (2003). PSORT-B: 
# improving protein subcellular localization prediction for 
# Gram-negative bacteria. Nuc Acids Res 31(13):3613-17. Please 
# cite this publication if you use PSORT-B in your research.

# The standalone version of PSORT-B is distributed under the GNU 
# General Public Licence (Gnu GPL) (see the LICENSE file included 
# in the download) by the Brinkman Laboratory, Simon Fraser 
# University, Burnaby, B.C., Canada.

# This standalone version of PSORT-B has initially been developed 
# for the Linux environment.

# This document describes the installation of the PSORT-B version 
# 1.1.4 command line program and the PSORT-B server packages. For 
# most purposes, following the installation instructions for the 
# command line version will be sufficient.

# For further information, please contact psort-mail@sfu.ca.


use Bio::Tools::Signal::Report;
use Bio::Tools::PSort::SVMLoc::DataSet;
use Bio::Tools::PSort::SVMLoc;
use Algorithm::HMM;

use Bio::Root::Root;

use strict;

# SubLoc localizations.
my %LOCS  = (Cytoplasmic   => 0, InnerMembrane => 1, Periplasmic   => 2,
	     OuterMembrane => 3, Extracellular => 4, Mitochondrial => 5,
	     Nuclear       => 6, Unknown       => 7);

my %LOCSR = (0 => 'Cytoplasmic',   1 => 'InnerMembrane', 2 => 'Periplasmic',
	     3 => 'OuterMembrane', 4 => 'Extracellular', 5 => 'Mitochondrial',
	     6 => 'Nuclear',       7 => 'Unknown');


use vars qw(@ISA $VERSION);
@ISA = qw(Bio::Root::Root);
$VERSION = '0.01';

# List of all the amino acid symbols.
my @AA = qw(V L I M F W Y G A P S T C H R K Q E N D);

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  # Get the model filenames and cutoff values, throwing an error if they
  # weren't provided.
  my ($svm, $hmm, $svmcutoff, $hmmcutoff)
    = $self->_rearrange(["SVM", "HMM", "SVMCUTOFF", "HMMCUTOFF"], @args);
  $self->throw("SVM model not specified") unless $svm;
  $self->throw("HMM model not specified") unless $hmm;

  # Create the SVM, throwing an error if we failed instantiate it.
  eval { $self->{svm} = new Bio::Tools::PSort::SVMLoc(Model => $svm); };
  $self->throw($@) if($@);

  # Create the HMM, throwing an error if we failed instantiate it.
  eval { $self->{hmm} = new Algorithm::HMM(Model => $hmm); };
  $self->throw($@) if($@);

  # Set the default cutoff value for the HMM and SVM.
  $self->{svmcutoff} = (defined($svmcutoff)) ? ($svmcutoff + 0) : 0.05;
  $self->{hmmcutoff} = (defined($hmmcutoff)) ? ($hmmcutoff + 0) : 0.01;

  return $self;
}

sub analyze {
  my ($self, $seq) = @_;
  my (%count, $subseq, $len, $ds, $loc, $res);

  # Make sure that we received a sequence object.
  $self->throw("Not a Bio::Seq object")
    if((! ref($seq)) || (! $seq->isa("Bio::SeqI")));

  $res = new Bio::Tools::Signal::Report();

  # We want at most the first 70 amino acids to analyze.
  $len = ($seq->length() > 70) ? 70 : $seq->length();
  $subseq = $seq->subseq(1, $len);

  # Run the sub-sequence through the HMM.
  my $rep = $self->{hmm}->search($subseq);
  my $hit = ($rep->domain_hits())[0];

  if($hit) {
    # Get the pvalue for the hit, and return immediately if we satisfy the
    # HMM cutoff.  (ie. We're pretty damn sure this looks like a cleavage
    # site.
    my $pval = $hit->pvalue;

    $res->hmm_score($hit->score);
    $res->hmm_pvalue($pval);
  }

  # Check for the lipoprotein motif.
  my $lip = $self->_is_prokar_lipoprotein($subseq);
  $loc = $lip ? $self->_localization($subseq, $lip) : 'Unknown';

  $res->is_lipoprotein($lip);
  $res->localization($loc);

  # Do an amino acid composition from the start of the sequence to five
  # amino acids after the predicted clevage point.
  my $ssubseq = substr($subseq, 0, $hit->seq_to);
  my $ssublen = length($ssubseq);

  # Count the occurences of each amino acid in the sub-sub-sequence.
  $count{uc($_)}++ for(split('', $ssubseq));

  # Create the dataset with all of the elements.
  $ds = new Bio::Tools::PSort::SVMLoc::DataSet(Label => 1);
  $ds->attribute($_, $count{$AA[$_ - 1]}/$ssublen || 0) for(1..scalar(@AA));

  $res->svm_prediction($self->{svm}->predict($ds));

  return $res;
}

sub _is_prokar_lipoprotein {
  my ($self, $seq) = @_;

  # Check for the motif (Prosite ID PS00013).
  if($seq =~ /[^.+DERK]{6}[LIVMFWSTAG]{2}[LIVMFYSTAGCQ][AGS]C/ig) {
    # Check to see where the cysteine residue was and whether or not
    # there was a charged residue (Lys or Arg) in the first seven residues.
    my $p = pos($seq);

    if((($p > 15) && ($p <= 36)) && (substr($seq, 0, 7) =~ /[RK]/)) {
      return $p if(substr($seq, $p + 1, 2) =~ /[DE]/);
    }
  }

  return 0;
}

sub _localization {
  my ($self, $seq, $pos) = @_;

  return substr($seq, $pos+1, 2) =~ /[DE]/ ? 'CytoplasmicMembrane': 'OuterMembrane';
}

sub train {
  my ($self, $data, $opts) = @_;
  my @tset;

  # Ensure we have at least 2 localizations to sort to.
  $self->throw("Training set must have > 1 localization")
    if(keys(%{$data}) < 2);

  # Create a new SVM.
  my $svm = new Bio::Tools::PSort::SVMLoc(Type => 'C-SVC', Kernel => 'linear');
  $self->throw("Error initializing SVM") if(! $svm);

  foreach my $loc (keys(%{$data})) {
    # Ensure we have a valid localization.
    $self->throw("Unknown localization: $_") if(! exists($data->{$loc}));

    # Ensure we have an array ref of something...
    $self->throw("Training set must be an array of Bio::Seq objects")
      if(ref($data->{$loc}) ne "ARRAY");

    for(@{$data->{$loc}}) {
       my (%count, $len, $ds);

       # Ensure we recieved a set of Bio::Seq objects.
       $self->throw("Not a Bio::Seq object") if(ref($_) ne "Bio::Seq");

       # Count the occurences of each amino acid in the sequence.
       $count{uc($_)}++ for(split('', $_->seq()));
       $len = $_->length();

       if($len) {
	 # Create a new dataset and add it to the list.
	 $ds = new Bio::Tools::PSort::SVMLoc::DataSet(Label => $LOCS{$loc});
	 $ds->attribute($_, ($count{$AA[$_-1]} || 0)/$len) for(1..scalar(@AA));

	 push(@tset, $ds);
       }
     }
  }

  $svm->train(@tset);

  $svm->save($opts->{Filename});
}

=head1 NAME

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 CONSTRUCTOR

=head1 METHODS

=head1 AUTHOR

Cory Spencer <cspencer@sfu.ca>

=head1 Maintainer

Matthew Laird <lairdm@sfu.ca>

=head1 SEE ALSO

=head1 ACKNOWLEGEMENTS

=cut

1;
