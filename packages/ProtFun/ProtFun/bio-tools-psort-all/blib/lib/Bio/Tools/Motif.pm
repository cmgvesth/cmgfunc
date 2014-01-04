package Bio::Tools::Motif;

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


use Bio::Tools::Motif::Pattern;
use Bio::Tools::Motif::Match;

use Bio::Root::Root;

use Carp;

use strict;

our @ISA = qw(Bio::Root::Root);
our $VERSION = '0.01';

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  # Get the model filename.
  my ($db) = $self->_rearrange(["DATABASE"], @args);

  $self->{patterns} = [ ];

  # Load database if specified.
  $self->load($db) if(defined($db));

  return $self;
}

sub match {
  my ($self, $seq) = @_;
  my @matches = ( );

  # Make sure that we received a sequence object.
  $self->throw("Not a Bio::Seq object")
    if((! ref($seq)) || (! $seq->isa("Bio::SeqI")));

  for(@{$self->{patterns}}) {
    my $regexp = $_->regexp;
    if($seq->seq =~ /$regexp/i) {
      my $start = length($`);
      my $end   = length($`) + length($&);

      my $m = new Bio::Tools::Motif::Match(-motifid      => $_->motif_id,
					   -localization => $_->localization,
					   -comment      => $_->comment,
					   -start        => $start,
					   -end          => $end);
      push(@matches, $m);
    }
  }

  return @matches;
}

sub add_pattern {
  my ($self, $pattern) = @_;

  $self->throw("Not a Bio::Tools::Motif::Pattern object")
    if((! ref($pattern) && (! $pattern->isa("Bio::Tools::Motif::Pattern"))));

  push(@{$self->{patterns}}, $pattern);
}

sub save {
  my ($self, $filename) = @_;

  open(NEWDB, ">$filename") || $self->throw("open: $filename: $!");
  for(@{$self->{patterns}}) {
    print($_->motif_id . "\t" . $_->localization . "\t" . $_->regexp ."\t" .
	  $_->comment ."\n");
  }
  close(NEWDB) || $self->throw("close: $filename: $!");

  return;
}

sub load {
  my ($self, $filename) = @_;

  open(DB, $filename) || $self->throw("open: $filename: $!");
  my @data = <DB>;
  close(DB) || $self->throw("close: $filename: $!");

  my $i = 1;
  for(@data) {
    chomp($_);
    my($id, $loc, $regexp, $comment) = split(/\t/, $_);

    if(defined($id) && defined($loc) && defined($regexp)) {
      my $pattern = new Bio::Tools::Motif::Pattern(-motifid      => $id,
						   -localization => $loc,
						   -regexp       => $regexp,
						   -comment      => $comment);

      push(@{$self->{patterns}}, $pattern);
    } else {
      $self->warn("line $i: suspicious looking motif database entry");
    }
  }

  return 1;
}

=head1 NAME

Bio::Tools::Motif - Perl implementation of the Motif protein subcellular
localization method.

=head1 SYNOPSIS

  use Bio::Tools::Motif;

  # Load a previously trained model from a file.
  $motif = new Bio::Tools::Motif(-database => 'motif-db.txt');

  # Attempt to match on a Bio::Seq object.
  @matches = $motif->match($seq);
  print($seq->display_id . ": matched " . $_->motif_id . "\n") for(@matches);

  # Save the database to a file.
  $motif->save('motif-db.txt.new');

  # Load a database from a file.
  $motif->load('motif-db.txt.new');

  # Add a new Bio::Tools::Motif::Pattern.
  $motif->add_pattern($pattern);

=head1 DESCRIPTION

Bio::Tools::Motif uses a selection of motifs that have been identified
to be typical of proteins resident at a specific subcellular localization.
The module accepts a Bio::Seq object and attempts to match it against
a database, returning one or more Bio::Tools::Motif::Match objects with
the prediction information if successful.

=head1 CONSTRUCTOR

   $motif = new Bio::Tools::Motif(-database => 'motif-db.txt');

The Motif constructor accepts the name of an existing database file.

=head1 METHODS

   @matches = $motif->match($seq);

The match method accepts a Bio::Seq object as an argument and returns
an array of Bio::Tools::Motif::Match objects that matched the given
sequence.

  $motif->add_pattern($pattern)

The add_pattern method adds another Bio::Tools::Motif::Pattern method
to the list of patterns used to match against sequences.

  $motif->save('motif-db.txt.new');

Saves the currently loaded database to the specified file.  Returns true on
success, false on failure.

  $motif->load('motif-db.txt.new');

Loads an existing database from file.  Returns true on success, false on
failure.

The database is a file containing a series of tab delimited fields.  The
fields (in order) are: the motif ID, the localization site, the perl
regular expression used to match the motif and an optional comment field.

=head1 AUTHOR

Cory Spencer <cspencer@sfu.ca>

=head1 SEE ALSO

Bio::Tools::Motif::Pattern, Bio::Tools::Motif::Match

=head1 ACKNOWLEGEMENTS

Thanks go out to Fiona Brinkman, Jennifer Gardy and the other members of the
Simon Fraser University Brinkman laboratory.

=cut

1;
