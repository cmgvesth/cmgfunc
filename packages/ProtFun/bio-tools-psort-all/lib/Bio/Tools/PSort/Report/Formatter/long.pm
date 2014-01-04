package Bio::Tools::PSort::Report::Formatter::long;

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
use Bio::Tools::PSort;

use vars qw(@ISA);
use strict;

@ISA = qw(Bio::Root::Root);

sub new {
  my ($class, @args) = @_;
  my $self = bless({ }, $class);

  ($self->{cutoff}, $self->{gram}, $self->{psort}, $self->{amods}, $self->{mcutoff}, $self->{skiplocalization}) = 
      $self->_rearrange(['CUTOFF','GRAM','PSORT', 'AMODS', 'MCUTOFF', 'SKIPLOCALIZATIONS'], @args);
  $self->{cutoff} = 0.75 if(! $self->{cutoff});
  $self->{mcutoff} = 0.5 if(! $self->{mcutoff});

  return $self;
}

sub format {
  my ($self, @reps) = @_;
  my $output;

  # Localizations to display
  my @alllocs;
  if($self->{gram} eq 'gramneg') {
      @alllocs = qw(Cytoplasmic CytoplasmicMembrane Periplasmic OuterMembrane Extracellular);
  } elsif ($self->{gram} eq 'grampos') {
      @alllocs = qw(Cytoplasmic CytoplasmicMembrane Cellwall Extracellular);
  } elsif ($self->{gram} eq 'archaea') {
      @alllocs = qw(Cytoplasmic CytoplasmicMembrane Cellwall Extracellular);
  } else {
      return "Formatter::long: We should never be here, invalid localization";
  }

  if(@reps) {
    my @amods;
    # This is a nasty hack because when you use the tree branching of
    # a pathway the modules don't all get run, so how do you then list all
    # the modules in the long format?  You don't.  So you have to pass them
    # in as parameters.  For local psort we pass in a pointer to the PSORT
    # object and just ask for the modules in the analysis path, for remote
    # PSort, since we don't have access to the PSort object we have to just
    # pass in a list of the modules.  Ugly, I know.
    if($self->{psort}) {
	@amods = $self->{psort}->get_modules("analysis");
    } elsif($self->{amods}) {
	@amods = @{$self->{amods}};
    } else {
	@amods = $reps[0]->get_modules("analysis");
    }
    my @omods = $reps[0]->get_modules("output");

    # Create the report header.
    $output = "SeqID";
    $output .= "\t$_" . "_Localization\t$_" ."_Details" for(@amods);
    for(@omods) {
      if($_ =~ /Bayes/i) {
	for (@alllocs) {
	  $output .= "\t$_" . "_Score";
	}
	$output .= "\tFinal_Localization\tFinal_Localization_Details\tFinal_Score\tSecondary_Localization";
	$output .= "\tPSortb_Version" if($Bio::Tools::PSort::VERSION_str);

      }
    }
    $output .= "\n";

    # Add a new line for every sequence.
    for my $rep (@reps) {
	my @secondlocs = ();
	# Make sure we got the right kind of object.
	$self->throw("Not a Bio::Tools::PSort::Report")
	    if((! ref($rep)) && (ref($rep) ne "Bio::Tools::PSort::Report"));

	# Get the sequence we're about to format the report for.
	$output .= $rep->seq->display_id() . " " . $rep->seq->desc();

	# Output the Analysis details.
	for (@amods) {
	    my @res = $rep->get_result("analysis", $_);

	    # Problem, what happens if we have multiple localizations
	    # across multiple $res.... you get the following ugly
	    # piece of code that merges them all but takes
	    # out any duplicates.  You also have to worry about
	    # modules not returning a result because they didn't
	    # run, ie. when SCL-BlastE returns a result.
	    my @locs = ();
	    
	    for my $res (@res) {
		if(@res) {
		    for my $r ($res->localization) {
			if(! (grep /$r/, @locs)) {
			    push @locs, $r;
			}
		    }
		    for my $sloc ($res->secondarylocalization) {
			if(! (grep /$sloc/, @secondlocs)) {
			    push @secondlocs, $sloc;
			}
		    }
		} else {
		    push @locs, "Unknown";
		}
	    }

	    # Do some module specific touch-ups...
	    @locs = "Non-Cytoplasmic"
		if(($locs[0]) && ($locs[0] ne "Unknown") && ($_ =~ /Signal/i));
	    
	    my @dtl = ();
	    for my $res (@res) {
		if($res) {
		    push @dtl, join(',', $res->details);
		}
	    }

	    $output .= "\t" . join(',', @locs) . "\t" . join(',', @dtl);
	}

	# Output the Final Prediction stuff.
	for (@omods) {
	    if($_ =~ /Bayes/i) {
		my @res = $rep->get_result("output", $_);

		# We need to pull out the scores in a specific order (to match the
		# header) so dump them into a hash.
		my $scr;
		$scr->{$_->localization} = sprintf("%4.2f", ($_->score * 10))
		    for(@res);

		for (@alllocs) {
		    $output .= "\t$scr->{$_}";
		}

		# Get the Final prediction.
		my ($fscr) = @res;
		my $fscore = $fscr->score * 10;
		my $loc = ($fscr->score >= $self->{cutoff}) ? $fscr->localization : "Unknown";
		my $details = '';

		if($self->{skiplocalization}) {
		    if(grep {m|^$loc?$|} @{$self->{skiplocalization}}) {
			$loc = 'Unknown';
			$fscore = 10;
			$details = '(Predicted localization does not exist in this organism)';
		    }
		}

		$details .= " (This protein may have multiple localization sites.)" 
		    if(($fscr->score >= $self->{mcutoff}) && ($fscr->score < $self->{cutoff}));

		$output .= "\t$loc\t$details\t" . (sprintf("%4.2f", $fscore));
	    }
	}

	$output .= "\t";
	if(@secondlocs) {
	    $output .= join(',', @secondlocs);
	}
	
	$output .= "\t";
	if($Bio::Tools::PSort::VERSION_str) {
	    $output .= "$Bio::Tools::PSort::VERSION_str\t";
	}

	$output .= "\n";
    }
  }

  return $output;
}

1;
