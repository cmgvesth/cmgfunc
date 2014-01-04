package Bio::Tools::Run::SCLBlast;

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

use Bio::Tools::Run::SCLBlast::Report;
use Bio::Tools::Run::SCLBlast::Hit;

use Bio::Tools::Run::StandAloneBlast;
use Bio::Root::Root;

use strict;

#our @ISA = qw(Bio::Tools::Run::StandAloneBlast);
our @ISA = qw(Bio::Root::Root);
our $VERSION = '0.02';
our %ENV;

sub new {
  my ($class, %args) = @_;

  my $self = bless {
      exact => $args{-exact} || 0,
      -database => $args{-database} || '',
      -blastdir => $args{-blastdir} || '',
  }, $class;

  $self->{blast} = new Bio::Tools::Run::StandAloneBlast(database => $self->{-database},
						program  => 'blastp',
						_READMETHOD => 'Blast',
						'F' => 'F');
  $self->throw("Failed to initialize SCLBlast") if(! $self->{blast});
  $self->{blast}->cleanup;

  # Set the default e value.
  $self->{blast}->expect(10e-2);

  $self->{parser} = sub {
      my $line = shift;

      # Extract the localization/description strings.
#    $line =~ /^.*\[(.*)\]\s+(.*)$/;
#    $line =~ /^\s*\[(\S+)\]\s+(.*)$/;
      $line =~ /^\s*\[(\S+)\]\s*(\[\S+\])?\s+(.*)$/;
      my ($loc, $sec, $dsc) = ($1, $2, $3);

      # Check to see if we have multiple localizations.
      $loc = [split(/\s+|\,/, $loc)] if($loc =~ /\,/);

      # Check to see if we have secondary localizations
      my $secloc = [];
      if($sec) {
	  $sec =~ /\[(\S+)\]/;
	  $secloc = $1;
	  $secloc = [split(/\s+|\,/, $secloc)] if($secloc =~ /\,/);
      }

    return ($loc, $secloc, $dsc);
  };

  return $self;
}

sub parser {
  my ($self, $parser) = @_;

  $self->throw("Not a CODE reference") if(ref($parser) ne "CODE");
  $self->{parser} = $parser;
}

sub blast {
  my ($self, $seq) = @_;
  my ($blast, $brep);
  my $qper;
  my $hper;
  my $fracID;
  my $expect;

  # Ensure we received a Bio::Seq object.
  $self->throw("Not a Bio::Seq object")
    if((! ref($seq)) || (! $seq->isa('Bio::Seq')));

  # Blast the sequence.
  $brep = $self->{blast}->blastall($seq);
  $self->{blast}->cleanup;
  return () if(! $brep);

  my $bres = $brep->next_result();

  # Create a list of the hits on the sequence.
  my @hits = ();

  while(my $hit = $bres->next_hit()) {
    my @hsps;
    while(my $hsp = $hit->next_hsp) {
      $qper = $hsp->length('hit')/$hsp->length('query') * 100;
      $expect = $hsp->evalue();

      # Check the evalue to make sure it's in range (stupid Blast bug)
      next if($expect > 1e-9);

      # If we're in exact match mode, skip those not within 1 basepairs
      # in similar length
      next if($self->{exact} && (abs($hit->length() - $seq->length()) > 1));

      # Add the HSP to the list if it meets our length criteria.
      push(@hsps, $hsp) if(($qper > 80) && ($qper < 120));
      $fracID = $hsp->frac_identical('total')*100;
    }

    if(@hsps) {
      my ($loc, $secloc, $desc) = $self->{parser}->($hit->description);

      if(!(($self->{exact}) && !( $qper == 100 &&  $fracID == 100))) {
        my $hit = new Bio::Tools::Run::SCLBlast::Hit(
	            '-localization' => $loc,
		    '-secondloc' => $secloc,
		    '-description' => ($self->{exact}?'100% ':'') . $hit->name . ": " . $desc,
		    '-algorithm' => 'blastp',
		    '-name' => 'sclblast',
		    '-length' => $hit->length,
		    '-hsps' => \@hsps);
        push(@hits, $hit);
      }
    }
  }


      return new Bio::Tools::Run::SCLBlast::Report(-hits => \@hits);

}

1;
