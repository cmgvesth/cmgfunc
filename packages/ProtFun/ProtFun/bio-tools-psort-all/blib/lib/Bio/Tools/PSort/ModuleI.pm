=head1 NAME

Bio::Tools::PSort::ModuleI - an interface class that all PSort modules should derive
implement.

=head1 SYNOPSIS

    package Bio::Tools::PSort::Module::Foo;

    use Bio::Tools::PSort::Module;

    use vars qw(@ISA);
    @ISA = qw(Bio::Tools::PSort::Module;

    ...
    ...

    1;

=head1 DESCRIPTION

The module object contains functionality common to all PSort modules.  This
class is not intended to be instantiated on it's own, it only provides 
functionality common to all PSort modules.


=head1 METHODS

    $result = $module->run($seq);

The run method must be overridden in all derived classes.  run() should accept
a Bio::Tools::Seq object, and perform whatever analysis the module provides on the
sequence, returning a bio::PSort::Scort object upon completion.

=head1 SEE ALSO

Bio::Tools::PSort, Bio::Tools::PSort::Pathway, Bio::Tools::PSort::Score

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

=head1 LICENCE

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

package Bio::Tools::PSort::ModuleI;

use Bio::Root::Root;

use vars qw(@ISA);
@ISA = qw(Bio::Root::Root);

use strict;

sub module_name {
  my $self = shift;

  $self->throw_not_implemented;
}

sub module_description {
  my $self = shift;

  $self->throw_not_implemented();
}

sub module_type  {
  my $self = shift;

 SWITCH: {
    $self->isa("Bio::Tools::PSort::Module::InputI") && return 'input';
    $self->isa("Bio::Tools::PSort::Module::AnalysisI") && return 'analysis';
    $self->isa("Bio::Tools::PSort::Module::OutputI") && return 'output';

    return 'unknown';
  }
}

sub run {
  my $self = shift;

  $self->throw_not_implemented();
}

sub _run {
  my ($self, $rpt) = @_;
  my $res;

  $res = $self->run($rpt);

  return $res;
}

1;
