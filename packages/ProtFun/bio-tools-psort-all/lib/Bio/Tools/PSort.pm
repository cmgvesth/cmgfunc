=head1 NAME

Bio::Tools::PSort - Module implementing a method for predicting protein
localization sites.

=head1 SYNOPSIS

    $psort = new Bio::Tools::PSort();

    $psort->install('Motif', 'Motif', { -database => $motifdb });

    $psort->classify($seq);

=head1 DESCRIPTION

The PSort module implements a method for predicting protein localization sites
as dicussed in the paper 'Expert System for Predicting the Protein
Localization Sites in Gram-Negative Bacteria' by Kenta Nakai and Minoru
Kanehisa.

=head1 CONSTRUCTOR

    $psort = new Bio::Tools::PSort();

Creates a new Bio::Tools::PSort object.

=head1 METHODS

    $psort->install('Motif', 'Motif', { -database => $motifdb });

The install method can be used to load a module after the PSort class has
been instantiated.  It accepts the name of the module to load, an alias used to
refer to the module by the program, and an optional hashref of parameters that
should be passed to the modules constructor.  The alias allows the same module to
be loaded more than once, but with different configuration options.

    $result = $psort->classify($seq);

The clasif method attempts to determine the likely protein localization for the
provided sequence.  A Bio::Tools::PSort::Report object is returned.


=head1 SEE ALSO

Bio::Tools::PSort::ModuleI, Bio::Tools::PSort::Score, Bio::Tools::PSort::Report,
Bio::Tools::PSort::Report::Formatter

=head1 AUTHOR

 Fiona Brinkman, Cory Spencer <psort-mail@sfu.ca>
 Brinkman Laboratory, Simon Fraser University, BC, Canada

=head1 MAINTAINER

 Matthew Laird <lairdm@sfu.ca>
 Brinkman Laboratory, Simon Fraser University, BC, Canada

=head1 OVERVIEW

PSORT-B is described in Gardy, J.L. et al (2003). PSORT-B: 
improving protein subcellular localization prediction for 
Gram-negative bacteria. Nuc Acids Res 31(13):3613-17. Please 
cite this publication if you use PSORT-B in your research.

The standalone version of PSORT-B is distributed under the GNU 
General Public Licence (Gnu GPL) (see the LICENSE file included 
in the download) by the Brinkman Laboratory, Simon Fraser 
University, Burnaby, B.C., Canada.

This standalone version of PSORT-B has initially been developed 
for the Linux environment.

This document describes the installation of the PSORT-B version 
1.1.4 command line program and the PSORT-B server packages. For 
most purposes, following the installation instructions for the 
command line version will be sufficient.

For further information, please contact psort-mail@sfu.ca.

=head1 LICENSE

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

=cut

package Bio::Tools::PSort;

use Bio::Tools::PSort::Pathway;
use Bio::Tools::PSort::Report;

use Bio::Root::Root;

use vars qw(@ISA $VERSION);
use strict;

@ISA = qw(Bio::Root::Root);
$VERSION = 3.0;

our $VERSION_str = 'PSORTb version 3.0';

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($modules) = $self->_rearrange([qw(MODULES)], @args);

  # Attempt to install all of the modules specified in the constructor.
  $self->install($_) for(@{$modules});

  return $self;
}

sub classify {
  my ($self, $seq, $path) = @_;

  # Make sure that we received a sequence object.
  $self->throw("Not a Bio::Seq object") if(! $seq->isa("Bio::SeqI"));

  # Create a report to hold the results for the sequence.
  my $report = new Bio::Tools::PSort::Report('-seq' => $seq);

  # Traverse each of the input, analysis and output pathways.
  my %res;
  foreach my $type ('input', 'analysis', 'output') {
    my $name = $path->{$type} || 'all';
    my $pway;

    # Locate the appropriate pathway.
    if($name eq 'all') {
      # If the pathway had the special name 'all', then dynamically create
      # a pathway consisting of all modules of the type.
      my @mods = keys %{$self->{modules}->{$type}};

      if(@mods) {
	my $name = shift @mods;
	my $path = [$self->{modules}->{$type}->{$name}, $name];

	$path = [$self->{modules}->{$type}->{$_}, $_, sub { 1 }, $path]
	  for (@mods);

	$pway = new Bio::Tools::PSort::Pathway('-path' => $path);
      }
    } else {
      $pway = $self->{pathways}->{$type}->{$name};
      $self->throw("$type pathway not found: $name") if(! defined($pway));
    }

    # Run the sequence through the pathway, passing in the results from the
    # previous pathway.
    %res = (defined($pway)) ? $pway->traverse($seq, %res) : ( );

    $report->add_results($type, %res);
  }

  # Return the final report.
  return $report;
}

sub install {
  my ($self, $mod, $alias, $args) = @_;

  # If we didn't receive the name or alias of a module to load, throw an
  # exception.
  $self->throw("No module specified") if(! $mod);
  $self->throw("No module alias specified") if(! $alias);

  $args = { } if(! defined($args));

  # Attempt to load the specified module.
  eval("require Bio::Tools::PSort::Module::$mod");
  $self->throw("Failed to load module $mod: $@") if($@);

  # Attempt to instantiate the module, throw exception on failure.
  my $module = eval { "Bio::Tools::PSort::Module::$mod"->new(%{$args}) };
  $self->throw("Failed to instantiate module $mod: $@") if($@);

  # Make sure the module we just created is actually a module.
  $self->throw("Not a Bio::Tools::PSort::ModuleI object")
    if(! $module->isa("Bio::Tools::PSort::ModuleI"));

  # Determine what type of module was loaded and place it in the appropriate
  # list(s).

  # Input module
  $self->{modules}->{input}->{$alias} = $module
    if $module->isa("Bio::Tools::PSort::Module::InputI");

  # Analysis module
  if($module->isa("Bio::Tools::PSort::Module::AnalysisI")) {
    $self->{modules}->{analysis}->{$alias} = $module;

    # Setup a default pathway that traverses only the module itself.
    my $pathway = new Bio::Tools::PSort::Pathway(-path => [$ module ]);
    $self->{pathways}->{analysis}->{$alias} = $pathway;
  }

  # Output module
  $self->{modules}->{output}->{$alias} = $module
    if $module->isa("Bio::Tools::PSort::Module::OutputI");

  # Throw an exception if we weren't derived from the three main types of
  # modules.
  $self->throw("Not an analysis, input or output module")
    if(! $module->isa("Bio::Tools::PSort::ModuleI"));

  return 1;
}

sub add_path {
  my ($self, $name, $type, $path) = @_;

  # Ensure we have a valid pathway type.
  $self->throw("invalid pathway type: $type")
    if(($type ne "input") && ($type ne "output") && ($type ne "analysis"));

  # Ensure that the path was well formed.
  $self->_prepare_path($path, $type);

  # Add the pathway.
   $self->{pathways}->{$type}->{$name} = new Bio::Tools::PSort::Pathway(-path => $path);

  return 1;
}

sub _prepare_path {
  my ($self, $path, $type) = @_;

  # Check to ensure the path is in the proper format.
  $self->throw("not an array reference")
    if(ref($path) ne "ARRAY");

  my ($node, $func, $right, $left) = @{$path};

  # Bail out if the module has not been loaded.
  $self->throw("module not installed: $node")
    if(! exists($self->{modules}->{$type}->{$node}));

  # Modify the path, putting the correct module onto the front of the array.
  unshift(@{$path}, $self->{modules}->{$type}->{$node});

  # Ensure we have an anonymous subroutine.
  $self->throw("not a code reference")
    if($func && (ref($func) ne "CODE"));

  # Check the right and left sides of the tree.
  $self->_prepare_path($right, $type) if $right;
  $self->_prepare_path($left, $type) if $left;


  return 1;
}

sub list_pathways {
  my $self = shift;
  my @pathways;

  foreach my $type ('input', 'analysis', 'output') {
    foreach my $path (keys(%{$self->{pathways}->{$type}})) {
      print("type = $type\n");
      print("name = $path\n");
    }
  }

  return @pathways;
}

sub list_modules {
  my $self = shift;
  my @modules;

  foreach my $type ('input', 'analysis', 'output') {
    foreach my $mod (values(%{$self->{modules}->{$type}})) {
      print("MOD = $mod\n");
      my @row = ($type, $mod->module_name, $mod->module_description);
      push(@modules, \@row);
    }
  }

  return @modules;
}

sub get_pathway {
    my ($self, $type) = @_;

    return sort keys(%{$self->{pathways}->{$type}});
}

sub get_modules {
    my ($self, $type) = @_;

    return sort keys(%{$self->{modules}->{$type}});
}

1;
