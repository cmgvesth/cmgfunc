package Bio::Tools::PSort::Report::Formatter::normal;

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

  ($self->{cutoff}, $self->{gram}, $self->{mcutoff}, $self->{skiplocalization}) = $self->_rearrange(['CUTOFF', 'GRAM', 'MCUTOFF', 'SKIPLOCALIZATIONS'], @args);
  $self->{cutoff} = 0.75 if(! $self->{cutoff});
  $self->{mcutoff} = 0.5 if(! $self->{mcutoff});

  return $self;
}

sub format {
  my ($self, @reps) = @_;
  my $output;

  if(@reps) {
    for my $rep (@reps) {
      my $seq = $rep->seq();
      my $ismult;
      my $hasSecloc = 0;
      my @secondlocs = ();

      $output .= "SeqID: " . $seq->display_id . " " . $seq->desc . "\n";
      $output .= "  Analysis Report:\n";

      for my $mod ($rep->get_modules("analysis")) {
	my @res = $rep->get_result("analysis", $mod);
	for my $res (@res) {
	  my $details = join(', ', $res->details) || "No details";
	  my @locs = $res->localization;

	  if($res->secondarylocalization) {
	      $hasSecloc = 1;
	      for my $sloc ($res->secondarylocalization) {
		  if(! (grep /$sloc/, @secondlocs)) {
		      push @secondlocs, $sloc;
			}
	      }
	  }

	  # Flag the result if it has multiple localization sites.
	  $ismult = 1 if((@locs > 1) && ($mod !~ /Signal/i));

	  # Totally nastly hack, but Jenn wanted it so I'm blaming her.
	  @locs = "Non-Cytoplasmic"
	    if(($locs[0] ne "Unknown") && ($mod =~ /Signal/i));
	  $output .= sprintf("    %-18s%-30s[%s]\n", $mod, join(', ', @locs),
			     $details);
	}
      }

      $output .= "  Localization Scores:\n";
      my @res = $rep->get_result("output", "Bayesian");
      for my $res (@res) {
	$output .= sprintf("    %-22s %-.2f\n", $res->localization, ($res->score * 10));
      }

      $output .= "  Final Prediction:\n";
      my $topscore = shift(@res);

      if($topscore->score >= $self->{cutoff}) {
	my $loc = $topscore->localization;
	my $scr = $topscore->score * 10;

	if($self->{skiplocalization}) {
	    if(grep {m|^$loc?$|} @{$self->{skiplocalization}}) {
		$loc = 'Unknown (Predicted localization does not exist in this organism)';
		$scr = 10;
	    }
	}

	$loc .= " (This protein may have multiple localization sites.)" if($ismult);
	$output .= sprintf("    %-22s %-.2f\n", $loc, $scr);
      } elsif($topscore->score >= $self->{mcutoff}) {
	$output .= "    Unknown (This protein may have multiple localization sites.)\n";
      } else {
	$output .= "    Unknown\n";
      }

      if(@secondlocs > 0) {
	  $output .= "  Secondary localization(s):\n";
	  $output .= sprintf("    %-22s\n", join(',', @secondlocs));
      }

      $output .= "\n" . ("-" x 79) . "\n\n";
    }
  }

  return $output;
}

1;
