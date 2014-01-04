package Bio::Tools::PSort::Report::Formatter::terse;

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

  ($self->{cutoff}, $self->{skiplocalization}) = $self->_rearrange(['CUTOFF', 'SKIPLOCALIZATIONS'], @args);
  $self->{cutoff} = 0.75 if(! $self->{cutoff});

  return $self;
}

sub format {
  my ($self, @reps) = @_;
  my $output;

  if(@reps) {
    my @mods = $reps[0]->get_modules("output");

    # Create the report header.
    $output = "SeqID";
#    $output .= "\t$_" . "_Localization\t$_" . "_Score" for(@mods);
    $output .= "\t" . "Localization\t" . "Score" for(@mods);
    $output .= "\n";

    # Add a new line for every sequence.
    for my $rep (@reps) {
      # Make sure we got the right kind of object.
      $self->throw("Not a Bio::Tools::PSort::Report")
	if((! ref($rep)) && (ref($rep) ne "Bio::Tools::PSort::Report"));

      # Get the sequence we're about to format the report for.
      $output .= $rep->seq->display_id() . " " . $rep->seq->desc();

      for (@mods) {
	my ($res) = $rep->get_result("output", $_);

	my $loc = ($res->score >= $self->{cutoff}) ? $res->localization : "Unknown";
	my $scr = $res->score * 10;

	if($self->{skiplocalization}) {
	    if(grep {m|^$loc?$|} @{$self->{skiplocalization}}) {
		$loc = 'Unknown (Predicted localization does not exist in this organism)';
		$scr = 10;
	    }
	}

	$output .= "\t$loc\t" . sprintf("%4.2f", $scr);
      }

      $output .= "\n";
    }
  }

  return $output;
}

1;
