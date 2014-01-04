package Bio::Tools::PSort::Module::Rules;

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


use Bio::Tools::PSort::Module::OutputI;
use Bio::Tools::PSort::Report::Result;
use Data::Dumper;

use vars qw(@ISA);
@ISA = qw(Bio::Tools::PSort::Module::OutputI);

use strict;

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  return $self;
}

sub run {
  my ($self, $seq, %res) = @_;
  my ($hassig);

  # Return HMMTOP prediction if one was made.
  if(defined($res{HMMTOP}) && (my ($hmmtop) = @{$res{HMMTOP}})) {
    $hassig = 1 if(join('', $hmmtop->details) =~ /signal peptide/i);
    return $hmmtop if($hmmtop->localization ne "Unknown");
  }

  # Return Motif prediction if one was made.
  if(defined($res{Motif}) && (my ($motif) = @{$res{Motif}})) {
    return $motif if($motif->localization ne "Unknown");
  }

  # Return OMPMotif prediction if one was made.
  if(defined($res{OMPMotif}) && (my ($omotif) = @{$res{OMPMotif}})) {
    return $omotif if($omotif->localization ne "Unknown");
  }

  # Return SCLBlast prediction if one was made.
  if(defined($res{SCLBlast}) && (my ($scl) = @{$res{SCLBlast}})) {
    return $scl if($scl->localization ne "Unknown");
  }

  # Check to see if the Signal module found a signal peptide.
  if(defined($res{Signal}) && (my ($signal) = @{$res{Signal}})) {
    $hassig = 1 if(join('', $signal->details) =~ /signal peptide/i);
  }

  # Return SubLoc prediction only if Cytoplasmic was predicted.
  if(defined($res{SubLoc}) && (my ($subloc) = @{$res{SubLoc}})) {
    return $subloc if(($subloc->localization eq "Cytoplasmic") && ! $hassig);
  }

  return new Bio::Tools::PSort::Report::Result(-loc => 'Unknown');
}

1;
