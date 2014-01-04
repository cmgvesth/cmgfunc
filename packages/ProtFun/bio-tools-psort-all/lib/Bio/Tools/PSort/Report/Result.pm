package Bio::Tools::PSort::Report::Result;

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

use vars qw(@ISA);
use strict;

@ISA = qw(Bio::Root::Root);

sub new {
  my ($class, @args) = @_;
  my $self = bless({ }, $class);

  my ($loc, $secondloc, $score, $details) =
    $self->_rearrange([qw(LOC SECONDLOC SCORE DETAILS)], @args);

  if(defined($loc)) {
    $self->{localization} = (ref($loc) eq "ARRAY") ? $loc : [ $loc ];
  } else {
    $self->{localization} = [ ];
  }

  if(defined($secondloc)) {
    $self->{secondarylocalization} = (ref($secondloc) eq "ARRAY") ? $secondloc : [ $secondloc ];
  } 
  else {
    $self->{secondarylocalization} = [ ];
  }

  $self->{details} = (defined($details)) ? (ref($details) eq "ARRAY") ? $details : [$details] : [ ];
  $self->{score} = (defined($score)) ? $score : 0;

  return $self;
}

sub localization {
  my ($self, $loc) = @_;

  if(defined($loc)) {
    $self->{localization} = (ref($loc) eq "ARRAY") ? $loc : [ $loc ];
  } else {
    return wantarray ? @{$self->{localization}} : @{$self->{localization}}[0];  }
}

sub secondarylocalization {
  my ($self, $secondloc) = @_;

  if(defined($secondloc)) {
    $self->{secondarylocalization} = (ref($secondloc) eq "ARRAY") ? $secondloc : [ $secondloc ];
  } 
  else {
    return wantarray ? @{$self->{secondarylocalization}} : @{$self->{secondarylocalization}}[0];  
  }
}

sub score {
  my $self = shift;

  return (@_) ? $self->{score} = shift : $self->{score};
}

sub details {
  my $self = shift;

  if(@_) {
    my $details = shift;
    $self->throw("Details must be an arrayref") if(ref($details) ne "ARRAY");
    $self->{details} = $details;
  } else {
    @{$self->{details}};
  }
}


1;
