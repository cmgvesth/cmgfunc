package Bio::Tools::PSort::Module::Null;

use Bio::Tools::PSort::Module::AnalysisI;
use Bio::Tools::PSort::Report::Result;

use vars qw(@ISA);
@ISA = qw(Bio::Tools::PSort::Module::AnalysisI);

use strict;

# Dummy module to do nothing.  Always returns Unknown

sub new {
    my($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    return $self;
}

sub run {
    my ($self, $seq) = @_;

    return new Bio::Tools::PSort::Report::Result(-loc => 'Unknown',
					  -details => ['Module skipped']);

}

1;

=head1 NAME

Null module

=head1 SYNOPSIS

This module is a place holder that will always return Unknown.

=head1 DESCRIPTION

=head1 CONSTRUCTOR

=head1 METHODS

=head1 SEE ALSO

=head1 AUTHOR
 
Fiona Brinkman, Matthew Laird <psort-mail@sfu.ca>
Brinkman Lab, Simon Fraser University, Canada

=cut
