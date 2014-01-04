package Bio::Tools::PSort::Module::HMMTOP;

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

use Bio::Tools::Run::HMMTOP;

use vars qw(@ISA);
@ISA = qw(Bio::Tools::PSort::Module::AnalysisI);

# The hmmtop processes have a habit of going zombie without this.
$SIG{CHLD} = 'IGNORE';

use strict;

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($dir, $loc, $cutoff) = $self->_rearrange([ qw(HMMTOPDIR LOC CUTOFF) ], @args);

  $self->{hmmtop} = new Bio::Tools::Run::HMMTOP(-hmmtopdir => $dir) ||
    $self->throw("Couldn't create new Bio::Tools::HMMTOP object");

  $self->throw("No localization specified in HMMTOP\n") if(!$loc);
  $self->{loc} = $loc;

  # Set the number of helices we expect for an CytoplasmicMembrane protein.
  $self->{cutoff} = (defined($cutoff) ? ($cutoff + 0) : 3);

  return $self;
}

sub run {
  my ($self, $seq, %prev) = @_;

  # Ensure we received a Bio::Seq object.
  $self->throw("Not a Bio::Seq object")
    if((! ref($seq)) && (! $seq->isa("Bio::Seq")));

  my $rep = $self->{hmmtop}->analyze($seq);
  my $num = $rep->num_helices;

  # Check to see if we have data from the Signal module - if we do, we should
  # reduce the helix count IF a signal peptide was detected and a TMH was found
  # before amino acid 70.
  if(exists($prev{Signal})) {
    my ($res) = @{$prev{Signal}};

    if($res->details =~ /^Signal peptide detected/) {
      for(rep->num_helices) {
	$num-- if($_->end <= 70);
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
