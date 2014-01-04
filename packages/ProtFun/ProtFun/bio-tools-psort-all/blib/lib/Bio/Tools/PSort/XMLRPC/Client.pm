package Bio::Tools::PSort::XMLRPC::Client;

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


use RPC::XML::Client;
use RPC::XML;

use Bio::Tools::PSort::Report::Result;
use Bio::Tools::PSort::Report;
use Bio::Root::Root;

use strict;

use vars qw(@ISA);
@ISA = qw(Bio::Root::Root);

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  # Get any of the parameters passed to the constructor.
  my ($server,$method) = $self->_rearrange([qw(SERVER METHOD)], @args);

  # Set some default values if they weren't received in the constructor.
  $self->{server} = $server || 'http://localhost/RPC';
  $self->{method} = $method || 'psort.classify';

  # Create a new XML-RPC client.
  $self->{client} = RPC::XML::Client->new($self->{server});

  return $self;
}

sub classify {
  my ($self, $seq, $path) = @_;
  my ($req, $rep, $res);

  # Create the structure of the request we're going to send.
  $req = { seq  => $seq->seq,
	   path => { input    => $path->{input}    || 'all',
		     analysis => $path->{analysis} || 'all',
		     output   => $path->{output}   || 'all' }};

  $res = $self->{client}->send_request($self->{method}, $req);

  # Check to see if we received a reference to an object from the send_request.
  # (If we didn't, then there was a client side send error - check parameters
  # passed to the RPC::XML::Client constructor.)
  if(ref($res)) {
    if(! $res->is_fault) {
      my $res = $res->value;
      my $rep = new Bio::Tools::PSort::Report(-seq => $seq);
      for my $path (keys(%{$res})) {
	my %res = ( );
	for my $mod (keys(%{$res->{$path}})) {
	  my @res;
	  for(@{$res->{$path}->{$mod}}) {
	    my @locs = split(',', $_->{localization});
	    push(@res, new Bio::Tools::PSort::Report::Result(
		         -loc => \@locs,
			 -score => $_->{score} || 0,
                         -details => $_->{details}));
	  }

	  $res{$mod} = \@res;
	}
	$rep->add_results($path, %res);
      }

      return $rep;
    } else {
      # Fault was returned, throw an exception.
      $self->throw("XML-RPC fault: " . $res->string);
    }
  } else {
    $self->throw("XML-RPC client side send error");
  }
}

1;
