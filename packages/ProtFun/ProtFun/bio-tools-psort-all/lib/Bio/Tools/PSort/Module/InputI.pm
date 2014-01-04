=head1 NAME

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 METHODS

=head1 SEE ALSO

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

package Bio::PSort::Module::InputI;

use Bio::PSort::ModuleI;

use vars qw(@ISA);
@ISA = qw(Bio::PSort::ModuleI);

use strict;

sub run {
  my $self = shift;

  $self->throw_not_implemented();
}

1;
