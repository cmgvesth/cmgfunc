package Bio::Tools::Run::SCLBlast::Hit;

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


use Bio::Search::Hit::GenericHit;
use Bio::Root::Root;

use strict;

our @ISA = qw(Bio::Search::Hit::GenericHit);

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($loc, $secloc, $desc) = $self->_rearrange([qw(LOCALIZATION SECONDLOC DESCRIPTION)], @args);
  $loc = [$loc] if(ref($loc) ne "ARRAY");
  $secloc = [$secloc] if(ref($secloc) ne "ARRAY");

  $self->{_locs} = $loc || $self->throw("Must provide a localization site.");
  $self->{_seclocs} = $secloc || '';
  $self->{_desc} = $desc || '';

  return $self;
}

sub num_locs {
  my $self = shift;

  return scalar(@{$self->{_locs}});
}

sub num_seclocs {
  my $self = shift;

  return (@_) ? scalar(@{$self->{_seclocs}}) : 0;
}

sub description {
  my $self = shift;

  return (@_) ? $self->{_desc} = shift : $self->{_desc};
}

sub localization {
  my $self = shift;

  if(@_) {
    my $loc = shift;
    $self->{_locs} = (ref($loc) eq "ARRAY") ? $loc : [$loc];
  } else {
    use Data::Dumper;
    return (wantarray) ? @{$self->{_locs}} : $self->{_locs}->[0];
  }
}

sub secondarylocalization {
   my $self = shift;

   if (@_) {
      my $secondloc = shift;
      $self->{_seclocs} = (ref($secondloc) eq "ARRAY") ? $secondloc : [$secondloc];
   } 
   else {
     use Data::Dumper;
     return (wantarray) ? @{$self->{_seclocs}} : $self->{_seclocs}->[0];
   }

}

1;
