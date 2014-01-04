=head1 NAME

Bio::Tools::PSort::ModHMM - Wrapper for the Prodiv_TMHMM transmembrane helix predictor.

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 CONSTRUCTOR

=head1 SEE ALSO

=head1 AUTHOR

 Matthew Laird <lairdm@sfu.ca>
 Brinkman Laboratory, Simon Fraser University, BC, CANADA

=head1 OVERVIEW

PSORT-B is described in Gardy, J.L. et al (2005). PSORT-B
V.2.0: Expanded prediction of bacterial protein subcellular
localization and insights gained from comparative proteome
analysis. Bioinformatics. 21:617-623.  Please cite this
publication if you use PSORT-B in your research.

The standalone version of PSORT-B is distributed under the GNU 
General Public Licence (Gnu GPL) (see the LICENSE file included 
in the download) by the Brinkman Laboratory, Simon Fraser 
University, Burnaby, B.C., Canada.

This standalone version of PSORT-B has initially been developed 
for the Linux environment.

This document describes the installation of the PSORT-B version 
2.x command line program and the PSORT-B server packages. For 
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

package Bio::Tools::PSort::ModHMM;

use strict;
use 5.006;
use Carp;

require DynaLoader;
require Exporter;
use AutoLoader;

our @ISA = qw(Exporter DynaLoader);

our %EXPORT_TAGS = ( 'all' => [ qw( ) ] );
our @EXPORT = qw( findHelices getHelices );

our $VERSION = '0.01';

sub AUTOLOAD {
  my $constname;
  our $AUTOLOAD;
  ($constname = $AUTOLOAD) =~ s/.*:://;
  croak "& not defined" if $constname eq 'constant';
  my $val = constant($constname, @_ ? $_[0] : 0);
  if ($! != 0) {
    if ($! =~ /Invalid/ || $!{EINVAL}) {
      $AutoLoader::AUTOLOAD = $AUTOLOAD;
      goto &AutoLoader::AUTOLOAD;
    }
    else {
      croak "Your vendor has not defined Bio::Tools::PSort::ModHMM macro $constname";
    }
  }
  {
    no strict 'refs';
    # Fixed between 5.005_53 and 5.005_61
    if ($] >= 5.00561) {
      *$AUTOLOAD = sub () { $val };
    }
    else {
      *$AUTOLOAD = sub { $val };
    }
  }
  goto &$AUTOLOAD;
}

bootstrap Bio::Tools::PSort::ModHMM $VERSION;


sub findHelices {
    my ($seq, $hmmfile, $repfile, $path) = @_;

    # Replace BioPerl wildcards/invalid letters with the X wildcard
    # understood by TMHMM.
    $seq =~ s/[O\-\.\*\?~]/X/gi;

    return find_helices($seq, $hmmfile, $repfile, $path);
}

sub getHelices {
    my ($seq, $hmmfile, $repfile, $path) = @_;

    # Replace BioPerl wildcards/invalid letters with the X wildcard
    # understood by TMHMM.
    $seq =~ s/[O\-\.\*\?~]/X/gi;

    return get_helices($seq, $hmmfile, $repfile, $path);
}

1;

__END__
