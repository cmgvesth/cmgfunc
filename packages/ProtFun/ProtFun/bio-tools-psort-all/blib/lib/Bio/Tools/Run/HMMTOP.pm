=head1 NAME

Bio::Tools::Run::HMMTOP - Wrapper for the HMMTOP transmembrane helix predictor.

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 CONSTRUCTOR

=head1 SEE ALSO

=head1 AUTHOR

 Cory Spencer <cspencer@sfu.ca>
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

package Bio::Tools::Run::HMMTOP;

use Bio::Tools::Run::HMMTOP::Report;
use Bio::Tools::Run::HMMTOP::Helix;

use Bio::Root::Root;
use Bio::Root::IO;

our @ISA = qw(Bio::Root::Root);
our $VERSION = "0.01";

use strict;

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  ($self->{dir}) = $self->_rearrange([ qw(HMMTOPDIR) ], @args);

  # Make sure the HMMTOP program files exist.
  $self->throw("HMMTOP directory not found") if(! -e $self->{dir});
  $self->throw("HMMTOP program file not found")
    if (! -e Bio::Root::IO->catfile($self->{dir}, "hmmtop"));
  $self->throw("HMMTOP architecture file not found")
    if (! -e Bio::Root::IO->catfile($self->{dir}, "hmmtop.arch"));
  $self->throw("HMMTOP PSV file not found")
    if (! -e Bio::Root::IO->catfile($self->{dir}, "hmmtop.psv"));

  # Set some of the environment variables that HMMTOP needs to run.
  $ENV{HMMTOP_ARCH} = Bio::Root::IO->catfile($self->{dir}, "hmmtop.arch");
  $ENV{HMMTOP_PSV}  = Bio::Root::IO->catfile($self->{dir}, "hmmtop.psv");

  return $self;
}

sub analyze {
  my ($self, $seq) = @_;

  # Ensure we got a Bio::Seq object.
  $self->throw("Not a Bio::Seq object")
    if((! ref($seq)) || (! $seq->isa('Bio::Seq')));

  # Create an object we'll use to get temporary file names.
  my $io = new Bio::Root::IO();

  # Path to the HMMTOP executable.
  my $prog = Bio::Root::IO->catfile($self->{dir}, "hmmtop");

  # Create an temporary file to store the sequence in.
  my ($fh, $name) = $io->tempfile();
  print($fh ">" . $seq->display_id . "\n" . $seq->seq . "\n");
  close($fh);

  # Run HMMTOP and grab the output.
  my $tseq = $seq->seq;
  open(HMMTOP, "$prog -if=$name -sf=FAS -pi=spred -is=pseudo |");
  my $output = <HMMTOP>;
  close(HMMTOP);

  my @helices;
  if($output) {
    $output  =~ s!^.+(IN|OUT)\s+!!;
    my @data = split(/\s+/, $output);

    if(@data) {
      my $num = shift(@data);

      for(1..$num) {
	my ($start, $end) = (shift @data, shift @data);
	my $hlx = new Bio::Tools::Run::HMMTOP::Helix(-start  => $start,
						     -end    => $end);

	push(@helices, $hlx);
      }
    }
  }

  return new Bio::Tools::Run::HMMTOP::Report(-helices => \@helices);
}

1;

