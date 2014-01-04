package Bio::Tools::PSort::Module::ModHMM;

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

use Bio::Tools::PSort::ModHMM;

use vars qw(@ISA);
@ISA = qw(Bio::Tools::PSort::Module::AnalysisI);

use strict;

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($hmmfile, $repfile, $path, $loc, $cutoff) = $self->_rearrange([ qw(HMMFILE REPFILE PATH LOC CUTOFF) ], @args);

  $self->throw("No localization specified in ModHMM\n") if(!$loc);
  $self->{loc} = $loc;

  $self->throw("No HMM file specified in ModHMM\n") if(!$hmmfile);
  $self->{hmmfile} = $hmmfile;

  $self->throw("No replacement letter file specified in ModHMM\n") if(!$repfile);
  $self->{repfile} = $repfile;

  $self->throw("Path specified in ModHMM\n") if(!$path);
  $self->{path} = $path;

  # Set the number of helices we expect for an CytoplasmicMembrane protein.
  $self->{cutoff} = (defined($cutoff) ? ($cutoff + 0) : 3);

  return $self;
}

sub run {
  my ($self, $seq, %prev) = @_;

  # Ensure we received a Bio::Seq object.
  $self->throw("Not a Bio::Seq object")
    if((! ref($seq)) && (! $seq->isa("Bio::Seq")));

  my $rep = getHelices($seq->seq, $self->{hmmfile}, $self->{repfile}, $self->{path});
  my $num = $rep->count;

  $self->throw("Error openning HMM file in ModHMM\n") if($num == -10);

  $self->throw("Error openning replacement letter file in ModHMM\n") if($num == -11);

  $self->throw("Error with sequence in ModHMM\n") if($num == -12);

  $self->throw("Error with ModHMM path\n") if($num == -13);

  # Check to see if we have data from the Signal module - if we do, we should
  # reduce the helix count IF a signal peptide was detected and a TMH was found
  # before amino acid 70.
  if(exists($prev{Signal})) {
    my ($res) = @{$prev{Signal}};

    if($res->details =~ /^Signal peptide detected/) {
      for(my $i=1;$i <= $rep->count; $i++) {
	$num-- if($rep->end($i) <= 70);
      }
    }
  }

  my $details;
  SWITCH: {
      ($num == 0) && do {
	$details = "No internal helices found";
	last SWITCH;
      };

      ($num == 1) && do {
	$details = "1 internal helix found";
	last SWITCH;
      };

      $details = "$num internal helices found";
    };

  if($num >= $self->{cutoff}) {
    return new Bio::Tools::PSort::Report::Result(-details => [ $details ],
						 -score   => 0,
						 -loc     => $self->{loc});
  } else {
    return new Bio::Tools::PSort::Report::Result(-details => [ $details ],
						 -score   => 0,
						 -loc     => 'Unknown');
  }
}

1;
