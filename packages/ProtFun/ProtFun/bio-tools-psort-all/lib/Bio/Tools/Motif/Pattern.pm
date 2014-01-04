package Bio::Tools::Motif::Pattern;

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


use Bio::Root::Root;
use Carp;

use strict;

our @ISA = qw(Bio::Root::Root);
our $VERSION = '0.01';

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  # Get the parameters we've been passed.
  my ($id,$loc,$regexp,$comment) =
    $self->_rearrange([qw(MOTIFID LOCALIZATION REGEXP COMMENT)], @args);

  $self->{localization} = $loc;
  $self->{motif_id}     = $id;
  $self->{regexp}       = $regexp;
  $self->{comment}      = $comment;

  return $self;
}

sub motif_id {
  my $self = shift;

  return (@_) ? $self->{motif_id} = shift : $self->{motif_id};
}

sub localization {
  my $self = shift;

  return (@_) ? $self->{localization} = shift : $self->{localization};
}

sub regexp {
  my $self = shift;

  return (@_) ? $self->{regexp} = shift : $self->{regexp};
}

sub comment {
  my $self = shift;

  return (@_) ? $self->{comment} = shift : $self->{comment};
}

=head1 NAME

  Bio::Tools::Motif::Pattern - Represents a pattern used by the
  Bio::Tools::Motif module.

=head1 SYNOPSIS

  use Bio::Tools::Motif::Pattern;

  # Create a new object.
  $m = new Bio::Tools::Motif::Pattern(-motifid => 'P07110',
                                      -localization => 'OuterMembrane',
                                      -regexp => '(GG.G.D.*){4}',
                                      -comment => 'FIMBRIAL_USHER Pattern');

  # Accessor methods.
  $val = $m->motif_id();
  $val = $m->localization();
  $val = $m->regexp();
  $val = $m->comment();

=head1 AUTHOR

Cory Spencer <cspencer@sfu.ca>

=head1 SEE ALSO

  Bio::Tools::Motif

=head1 ACKNOWLEGEMENTS

Thanks go out to Fiona Brinkman, Jennifer Gardy and the other members of the
Simon Fraser University Brinkman laboratory.

=cut

1;
