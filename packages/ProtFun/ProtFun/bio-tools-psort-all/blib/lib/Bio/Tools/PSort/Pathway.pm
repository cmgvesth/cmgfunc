=head1 NAME

Bio::Tools::PSort::Pathway - an object used to define the pathway that is
traversed while doing a psort analysis.

=head1 SYNOPSIS

    $pathway = new Bio::Tools::PSort::Pathway('-path' => [ ... ]);

    @results = $pathway->traverse($seq);

=head1 DESCRIPTION

The Pathway object is used to describe the pathway that is traversed while
performing a psort analysis.  Traditional psort methods [1] involved
traversing a binary tree, where at every node an analysis was performed,
and based on the result, control was directed down one of the two branches.
This version of PSort maintains the ability to define binary tree pathways,
but allows them to be defined at runtime and/or loaded from a configuration
file.

Since the Pathway modules interface for defining new trees is rather low
level, it is recommended that end users use either of the
Bio::Tools::PSort::Pathway::Perl or Bio::Tools::PSort::Pathway::XML modules
to create a binary tree.

[1] As described in the paper 'Expert System for Predicting the Protein
Localization Sites in Gram-Negative Bacteria' by Kenta Nakai and Minoru
Kanehisa.

=head1 CONSTRUCTOR

    $pathway = new Bio::Tools::PSort::Pathway('-path' => [ ... ]);

The constructor for the Pathway object accepts a single named parameter,
'-path', which is an arrayref containing a description of the pathway.
Each node in the pathway is represented by an array ref containing
1 to 4 elements, which are in order, a reference to a
Bio::Tools::PSort::ModuleI subclass, a reference to a subroutine, and two
other node structures which will be executed if the subroutine returns either
true or false values.

Upon arriving at a node in the tree, the module is first run on the sequence
being analyzed, and it's results is passed to the subroutine.  If the
subroutine returns true, then the traversal continues with the third
parameter, another tree node, otherwise traversal is sent down the fourth
parameter, which is also another node in the tree.

For example:

    $mod1 = new Bio::Tools::PSort::Module::Regexp();
    $mod2 = new Bio::Tools::PSort::Module::SCLBLAST();
    $mod3 = new Bio::Tools::PSort::Module::SPaseI();

    $func = sub {
              my $res = shift;

              return ($res->score()->{cyoplasm} > 0.50);
            };

    $path = [$mod1, $func, [$mod2], [$mod3]];

    $pathway = new Bio::Tools::PSort::Pathway('-path' => $path);

This example describes a simple tree with three nodes.  The first node
consists of the Regular Expression module.  If the result returned from
the Regexp module has score higher that 0.50 for the cytoplasm localization,
then the subroutine will return true, and control flow will be directed
to the SCLBLAST module.  If the score is less than or equal to 0.50, the
control flow will be sent to the SPaseI module.

=head1 METHODS

    @results = $pathway->traverse($seq);

The traverse method is used to send a sequence down the pathway.  An
array of Bio::Tools::PSort::Score object is returned (each module returns a
single Score object, they are collected as the pathway is traversed
and returned as a group).

=head1 SEE ALSO

Bio::Tools::PSort, Bio::Tools::PSort::ModuleI, Bio::Tools::PSort::Score,
Bio::Tools::PSort::Pathway::Perl, Bio::Tools::PSort::Pathway::XML

=head1 AUTHOR

 Fiona Brinkman, Cory Spencer <psort-list@sfu.ca>
 Brinkman Laboratory, Simon Fraser University, BC, CANADA

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


package Bio::Tools::PSort::Pathway;

use Bio::Tools::PSort::Constants ':all';
use Bio::Root::Root;

use vars qw(@ISA);
use strict;

@ISA = qw(Bio::Root::Root);

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($path) = $self->_rearrange([qw(PATH)], @args);

  $self->{path} = $self->_build_path(@{$path});

  return $self;
}

sub traverse {
  my $self = shift;

  return $self->{path}->(@_);
}

sub _build_path {
  my $self = shift;
  my ($module, $name, $cond, $true, $false) = @_;

  # Make sure the module we're getting passed is actually a module.
  $self->throw("Not a Bio::Tools::PSort::ModuleI object")
    if(ref($module) && (! $module->isa("Bio::Tools::PSort::ModuleI")));

  # We got four parameters, build a complete node.
  if(defined($false)) {
      my $fun1 = $self->_build_path(@{$true});
      my $fun2 = $self->_build_path(@{$false});

      return sub($\%\%) {
	my @res = $module->run(@_);

	# Check the constant returned from the module - if it's OK, continue
	# processing the chain.  If it's not okay, return the constant up to
	# the top of the pathway.
	if($cond->(@res)) {
	  return ($name => \@res, $fun1->(@_, $name => \@res));
	} else {
	  return ($name => \@res, $fun2->(@_, $name => \@res));
	}
      };
    };

  # We got three parameters, branch down the left side of the tree if
  # $cond returns true, otherwise just return the result from the
  # current module.
  if(defined($true)) {
      my $fun = $self->_build_path(@{$true});

      return sub($\%\%) {
	my @res = $module->run(@_);

	if($cond->(@res)) {
	  return ($name => \@res, $fun->(@_, $name => \@res));
	} else {
	  return $name => \@res;
	}
      };
    };

  # We got either one or two parameters, just return the result from the
  # current module.
  if(defined($cond) || defined($module)) {
      return sub($\%\%) {
	my @res = $module->run(@_);

	return $name => \@res;
      };
    };
}

1;
