package Bio::Tools::PSort::Report::Formatter::html;

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

  # Localizations to display
  my @alllocs;
  if($self->{gram} eq 'gramneg') {
      @alllocs = qw(Cytoplasmic CytoplasmicMembrane Periplasmic OuterMembrane Extracellular);
  } else {
      @alllocs = qw(Cytoplasmic CytoplasmicMembrane Cellwall Extracellular);
  }

      $output .= "SeqID: " . $seq->display_id . " " . $seq->desc . "\n";
      $output .= "  Analysis Report:\n";

      for my $mod ($rep->get_modules("analysis")) {
	my @res = $rep->get_result("analysis", $mod);
	for my $res (@res) {
	  my $details = join(', ', $res->details) || "No details";
	  my @locs = $res->localization;

	  # Flag the result if it has multiple localization sites.
          $ismult = 1 if((@locs > 1) && ($mod !~ /Signal/i));

	  if($res->secondarylocalization) {
	      $hasSecloc = 1;
	      for my $sloc ($res->secondarylocalization) {
		  if(! (grep /$sloc/, @secondlocs)) {
		      push @secondlocs, $sloc;
			}
	      }
	  }

	  # Do some per-module output modifications.
	SWITCH:{
	    ($mod =~ /SCL\-?BLAST/i) && do {
	      if($details =~ /^(matched )(100% )?([OPQ]\d+):/) {
		$details =~ s/^(matched )(100% )?([OPQ]\d+):/$1$2<a href=\"http:\/\/www\.expasy\.org\/cgi\-bin\/niceprot\.pl\?$3\" target=\"blank\">$3<\/a>:/;
	      } elsif($details =~ /^(matched )(100% )?([\d.]+):/) {
		$details =~ s/^(matched )(100% )?([\d.]+):/$1$2<a href=\"http:\/\/www\.ncbi\.nlm\.nih\.gov\/entrez\/query\.fcgi\?cmd=Retrieve&db=protein&list_uids=$3&dopt=GenPept\" target=\"blank\">$3<\/a>:/;
	      }
	      last SWITCH;
	    };

	    ($mod =~ /^Motif/i) && do {
	      $details =~ s/^(matched )(PS\w+):/$1<a href=\"http:\/\/www\.expasy\.org\/cgi\-bin\/prosite\-search\-ac\?$2\" target=\"blank\">$2<\/a>/;
	      last SWITCH;
	    };

	    ($mod =~ /^Profile/i) && do {
	      $details =~ s/^(matched )(PS\w+):/$1<a href=\"http:\/\/www\.expasy\.org\/cgi\-bin\/prosite\-search\-ac\?$2\" target=\"blank\">$2<\/a>/;
	      last SWITCH;
	    };

	    ($mod =~ /^Signal/i) && do {
	      @locs = "Non-Cytoplasmic" if($locs[0] ne "Unknown");
	      last SWITCH;
	    }
	  }

	  $output .= sprintf("    %-18s%-30s[%s]\n", $mod, join(', ', @locs),
			     $details);
	}
      }

      $output .= "  Localization Scores:\n";
      my @res = $rep->get_result("output", "Bayesian");# ||
      my @res;
      if($rep->get_result("output", "Bayesian")) { @res = $rep->get_result("output", "Bayesian"); }
      elsif($rep->get_result("output", "Bayesian-")) { @res = $rep->get_result("output", "Bayesian-"); }
      elsif($rep->get_result("output", "Bayesian+")) { @res = $rep->get_result("output", "Bayesian+"); }
      elsif($rep->get_result("output", "Bayesian_a")) { @res = $rep->get_result("output", "Bayesian_a"); }
#	        $rep->get_result("output", "bayesn") ||
#	        $rep->get_result("output", "Bayesianp");
      my %results;
      # Ugly hack because Jenn and Sebastien wanted the the localizations
      # in a specific order, blame them.
      for my $res (@res) {
	  $results{$res->localization} = $res->score;
      }
      $output .= sprintf("    %-22s %-.2f\n", $_, ($results{$_} * 10))
	  for (@alllocs);

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

      $output .= "-" x 79 . "\n";
    }
  }

  return $output;
}

1;
