package Bio::Tools::PSort::Module::Profile;

use Bio::Tools::PSort::Module::AnalysisI;
use Bio::Tools::PSort::Report::Result;
use Bio::Tools::PSort::Profile;

use vars qw(@ISA);
@ISA = qw(Bio::Tools::PSort::Module::AnalysisI);

use strict;

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  #profile module must receive: 1 the path to the prosite profiles 2. the path to pfscan, 3. the path to the file containing the SCL's for each of the profile_ids
  my ($profiledb, $pfscan, $profileids) = $self->_rearrange([qw(DATABASE PROGRAM PROFILEIDS)], @args);

  #pass up the paths to Bio::Tools::Profile.
  $self->{profile} = new Bio::Tools::PSort::Profile(-database => $profiledb, 
                                             -program => $pfscan,
                                             -profileids => $profileids) ||
    $self->throw("Couldn't create new Bio::Tools::PSort::Profile object");

  return $self;
}


sub run {

  my ($self, $seq) = @_;

  my @res = map {
    my $detail = "matched " . $_->profile_id . ": " . $_->comment;
    new Bio::Tools::PSort::Report::Result(-details => [ $detail ],
					  -score   => 0,
					  -loc     => $_->localization);
  } $self->{profile}->match($seq);

  @res ? @res : new Bio::Tools::PSort::Report::Result(-loc => 'Unknown',
						      -details => ['No matches to profiles found']);
}

1;

=head1 NAME

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 CONSTRUCTOR

=head1 METHODS

=head1 SEE ALSO

=head1 AUTHOR
 
Fiona Brinkman, Chris Walsh, Matthew Laird <psort-mail@sfu.ca>
Brinkman Lab, Simon Fraser University, Canada

=cut

__END__
