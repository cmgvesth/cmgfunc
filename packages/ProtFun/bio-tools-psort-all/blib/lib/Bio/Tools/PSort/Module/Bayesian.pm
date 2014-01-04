package Bio::Tools::PSort::Module::Bayesian;

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

  my ($mfile, $prior) = $self->_rearrange([qw(MODEL PRIOR)], @args);

  $self->throw("Bayesian model file not found: $mfile") if(! -e $mfile);

  # Load the Bayesian model produced by the training program.
  open(MFILE, $mfile) || $self->throw("error: open $mfile: $!");
  my $model = join('', <MFILE>);
  close(MFILE) || $self->throw("error: close $model: $!");

  $self->{model} = eval($model);
  $self->throw("Error loading Bayesian model: $@") if($@);

  # Store any prior probabilities provided to us.  (If none provided,
  # each localization will be assumed to have an equal chance of
  # appearing.)
  if(defined($prior)) {
    $self->throw("Prior probabilties must be a hash reference")
      if(ref($prior) ne "HASH");
    $self->{prior} = $prior;
  } else {
    # If we weren't provided a list of prior probabilities, figure out all the
    # localizations we're going to be making predictions for.
    for my $anal (keys(%{$self->{model}})) {
      for my $loc (keys(%{$self->{model}->{$anal}})) {
	$self->{prior}->{$loc} = 1;
      }
    }
  }

  return $self;
}

sub run {
  my ($self, $seq, %res) = @_;
  my (%preds, @res);

  # Initialize to the prior probabilities of getting each localization.
  %preds = %{$self->{prior}};

  # This deals with SCL-BLASTe and 100% protein matches
  if($res{'SCL-BLASTe-'} || $res{'SCL-BLASTe+'} || $res{'SCL-BLASTe_a'}) {
      my $res;
      if($res{'SCL-BLASTe-'}) {
	  ($res) = @{$res{'SCL-BLASTe-'}};
      } elsif($res{'SCL-BLASTe+'}) {
	  ($res) = @{$res{'SCL-BLASTe+'}};
      } else {
	  ($res) = @{$res{'SCL-BLASTe_a'}};
      }
      my @locs = $res->localization;
      if(!($locs[0] eq 'Unknown')) {
	  my %scores = ();
	  $scores{$_} = 1
	      for(@locs);
	  push(@res, new Bio::Tools::PSort::Report::Result(-loc => $_, -score => ($scores{$_}?1:0)))
	      for(keys(%preds));

	  return sort({$b->score <=> $a->score} @res);
      }
  }

  for my $anal (keys(%res)) {
    my (@res) = @{$res{$anal}};

    # Ensure that the results from each analysis are the correct object type.
    for my $r (@res) {
      $self->throw("Not a Bio::Tools::PSort::Report::Result object")
	if((! ref($r)) && (! $r->isa("Bio::Tools::PSort::Report::Result")));
    }

    for my $res (@res) {
      # Get the prediction for the analysis.
      my @locs = $res->localization;

      for my $pred (@locs) {
	next if($pred eq "Unknown");
	# Update the probabilities for each location.
	$preds{$_} *= $self->{model}->{$anal}->{$_}->{$pred}
	  for(keys(%preds));
      }
    }
  }

  # Scale the scores so they all sum to 1.0
  my ($tot, $scale);
  $tot += $_ for(values(%preds));
  if($tot) {
    $scale = 1/$tot;
    $preds{$_} *= $scale for(keys(%preds));
  }

  # Sort the scores so the highest comes first.
  push(@res, new Bio::Tools::PSort::Report::Result(-loc => $_, -score => $preds{$_}))
    for(keys(%preds));

  return sort({$b->score <=> $a->score} @res);
}

1;

