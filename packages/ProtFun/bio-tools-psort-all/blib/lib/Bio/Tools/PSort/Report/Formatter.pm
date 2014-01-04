package Bio::Tools::PSort::Report::Formatter;

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

use vars qw/@ISA/;
use strict;

@ISA = qw(Bio::Root::Root);

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  return $self;
}

sub format {
  my ($self, $format, $args, @reps) = @_;

  $self->install($format, $args);

  return $self->{modules}->{$format}->format(@reps);
}

sub install {
  my ($self, $mod, $args) = @_;


  # If we didn't receive the name of a module to load, throw an exception.
  $self->throw("No module specified") if(! $mod);

  # If the module has already been loaded, return right away.
  return if(exists($self->{modules}->{$mod}));

  $args = { } if(! defined($args));

  # Attempt to load the specified module.
  eval("require Bio::Tools::PSort::Report::Formatter::$mod");
  $self->throw("Failed to load module $mod: $@") if($@);

  # Attempt to instantiate the module, throw exception on failure.
  my $module = eval { "Bio::Tools::PSort::Report::Formatter::$mod"->new(%{$args}) };
  $self->throw("Failed to instantiate module $mod: $@") if($@);

  # Store the module for later use.
  $self->{modules}->{$mod} = $module;

  return 1;
}

1;
