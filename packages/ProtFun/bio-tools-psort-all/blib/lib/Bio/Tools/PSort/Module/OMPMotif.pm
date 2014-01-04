package Bio::Tools::PSort::Module::OMPMotif;

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


use Bio::Tools::PSort::Module::AnalysisI;
use Bio::Tools::PSort::Report::Result;
use Bio::Tools::Motif;

use vars qw(@ISA);
@ISA = qw(Bio::Tools::PSort::Module::AnalysisI);

use strict;

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($db, $cutoff) = $self->_rearrange([qw(DATABASE CUTOFF)], @args);

  $self->{motif} = new Bio::Tools::Motif(-database => $db) ||
    $self->throw("Couldn't create new Bio::Tools::Motif object");

  # Set the number of patterns needed to make a classification, or default
  # to three if a value wasn't provided.
  $self->{cutoff} = defined($cutoff) ? ($cutoff + 0) : 3;

  return $self;
}

sub run {
  my ($self, $seq) = @_;

  my @matches = $self->{motif}->match($seq);
  if(@matches >= $self->{cutoff}) {
    my $detail = "matched " . @matches . " rules";
    $detail .= " (" . join(', ', map { $_->motif_id } @matches) . ")";

    return new Bio::Tools::PSort::Report::Result(-details => [ $detail ],
						 -score   => 0,
						 -loc     => 'OuterMembrane');
  }

  return new Bio::Tools::PSort::Report::Result(-loc => 'Unknown',
					       -details => ['No motifs found']);
}

1;

=head1 NAME

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 CONSTRUCTOR

=head1 METHODS

=head1 SEE ALSO

=head1 AUTHOR

 Cory Spencer <cspencer@sfu.ca>

=cut

__END__
