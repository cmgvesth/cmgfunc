package Bio::Tools::PSort::Module::Signal;

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


use Bio::Tools::PSort::Constants qw(:all);
use Bio::Tools::PSort::Module::AnalysisI;
use Bio::Tools::PSort::Report::Result;

use vars qw(@ISA);
use strict;

@ISA = qw(Bio::Tools::PSort::Module::AnalysisI);

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($svm, $hmm, $prog, $gram) = $self->_rearrange([ qw(SVM HMM PROGRAM GRAM) ], @args);

  $self->{svm} = $svm || $self->throw("No SVM model specified.");
  $self->{hmm} = $hmm || $self->throw("No HMM model specified.");
  $self->{program} = $prog || $self->throw("No signal program specified. ");
  $self->{gram} = ($gram eq 'Positive'?1:0);

  return $self;
}

sub run {
  my ($self, $seq, %prev) = @_;
  my ($hassig);

  $self->throw("Not a Bio::Seq object")
    if((! ref($seq)) && (! $seq->isa("Bio::Seq")));

  # Check to see if we can detect a signal peptide with the HMM/SVM module.
  my $sq = $seq->seq;

  # Run the program to check for signal peptides.
  open(SIG, "$self->{program} --svm=$self->{svm} --hmm=$self->{hmm} $sq |");
  my $res = <SIG>;
  chomp($res);
  close(SIG);

  $hassig = 1 if($res);

  # If a signal peptide was detected the sequence is one of the non-cytoplasmic
  # localizations.
  if($hassig) {
    my @locs;
    if($self->{gram}) {
	@locs = qw(CytoplasmicMembrane Extracellular Cellwall);
    } else {
	@locs = qw(CytoplasmicMembrane Periplasmic OuterMembrane Extracellular);
    }
    my $detail = "Signal peptide detected";
    return new Bio::Tools::PSort::Report::Result(-loc => \@locs,
						 -details => [ $detail ]);
  } else {
    my $detail = "No signal peptide detected";
    return new Bio::Tools::PSort::Report::Result(-loc => "Unknown",
						 -details => [ $detail ]);
  }
}

1;
