<?xml version="1.0"?>
<!DOCTYPE methoddef SYSTEM "rpc-method.dtd">
<!--
    Generated automatically by make_method v1.09, Mon Jul 14 09:48:57 2003

    Any changes made here will be lost.
-->
<methoddef>
<name>psort.classify</name>
<version>1.0</version>
<signature>struct struct</signature>
<help>
PSort routine.
</help>
<code language="perl">
<![CDATA[
#!/usr/bin/perl
use Apache::Bio::Tools::PSort;
use Bio::Seq;
use Data::Dumper;

use strict;

sub classify {
  my ($srv, $req) = @_;

  # Pull the sequence and path information out of the request.
  my $path = {input    => $req->{path}->{input}    || 'all',
	      analysis => $req->{path}->{analysis} || 'all',
	      output   => $req->{path}->{output}}  || 'all';
  my $seq  = new Bio::Seq(-seq => $req->{seq});

  # Get an instance of the PSort object.
  my $psort = Apache::Bio::Tools::PSort->instance();
  my $rep   = $psort->classify($seq, $path);

  my $res = { };
  for my $path qw(input analysis output) {
    my @mods = $rep->get_modules($path);

    for my $mod (@mods) {
      $res->{$path}->{$mod} = [ ];
      for($rep->get_result($path, $mod)) {
	my $r = {
		 localization => join(',', $_->localization),
		 score        => $_->score,
		 details      => [ $_->details ]
		};

	push(@{$res->{$path}->{$mod}}, $r);
      }
    }
  }

  return $res;
}

__END__
]]></code>
</methoddef>
