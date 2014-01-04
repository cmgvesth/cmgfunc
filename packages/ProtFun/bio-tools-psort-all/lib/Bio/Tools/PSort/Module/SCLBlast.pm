package Bio::Tools::PSort::Module::SCLBlast;

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


use Bio::Tools::PSort::Report::Result;
use Bio::Tools::PSort::Constants qw(:all);
use Bio::Tools::PSort::Module::AnalysisI;
use Bio::Tools::Run::SCLBlast;
use Bio::Tools::Run::SCLBlastLocal;
use Bio::SearchIO;

use vars qw(@ISA);
use strict;

@ISA = qw(Bio::Tools::PSort::Module::AnalysisI);

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($db, $exact, $blastbin) = $self->_rearrange([qw(DATABASE EXACT BLASTBIN)], @args);

  # Save the parameters for the reload option
  $self->{db} = $db;
  $self->{exact} = $exact;
  $self->{blastbin} = $blastbin;

  # Create a new Bio::Tools::Run::SBLBlast instance.
  if($blastbin) {
      $self->{sclblast} = new Bio::Tools::Run::SCLBlastLocal( -database => $db, -exact => $exact, -blastbin => $blastbin );
  } else {
      $self->{sclblast} = new Bio::Tools::Run::SCLBlast( -database => $db, -exact => $exact );
  }

  return $self;
}

sub run {
  my ($self, $seq) = @_;
  my @res;

  # Ensure that we received a Bio::Seq object.
  $self->throw("Didn't receive a Bio::Seq")
    if((! ref($seq)) || (! $seq->isa('Bio::Seq')));

  if(my $brep = $self->{sclblast}->blast($seq)) {
    my $hit = $brep->next_hit;
    if($hit) {
      my @locs = $hit->localization;
      my @seclocs = $hit->secondarylocalization;

      my $r = new Bio::Tools::PSort::Report::Result();
      $r->localization(\@locs);
      $r->secondarylocalization(\@seclocs);
      $r->details(["matched " . $hit->description]);
      push(@res, $r);
    }
  }

  @res ? @res : new Bio::Tools::PSort::Report::Result(-loc => 'Unknown',
						      -details => ['No matches against database']);
}

1;
