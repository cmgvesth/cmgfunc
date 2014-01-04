
package Bio::Tools::PSort::Profile;


use Bio::Tools::PSort::Profile::Match;
use Bio::SeqIO;
use Bio::Root::IO;
use Bio::Root::Root;

use Carp;

use strict;

our @ISA = qw(Bio::Root::Root);
our $VERSION = '0.01';

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);

    # Get the file paths for the profile and for the pfscan binary
    my ($profiledb, $pfscan, $profileids) = 
	$self->_rearrange([qw(DATABASE PROGRAM PROFILEIDS)], @args);

    $self->throw("No profile database selected") if (!$profiledb); 
    $self->{profile}->{profiledb} = $profiledb;
 
    $self->throw("No pfscan path selected") if (!$pfscan);
    $self->{profile}->{pfscan} = $pfscan;

    $self->throw("No profile ids file selected") if (!$profileids);
    #this file contains a the profile_ids and the corresponding localization
    open(PROFILEIDS, $profileids) or $self->throw("Could not open file $profileids");

    ###build a hash with ids as keys and scl's as values
    my %localization;

    while(<PROFILEIDS>){
	if (/^(\S+)\s+(\S+)/){
	    my $id = $1;
	    my $scl = $2; #'
	    chomp($scl);
	    $localization{$id} = $scl;
	}
    }
    close PROFILEIDS;
    $self->{profile}->{profileids} = \%localization;
 
    return $self;
}

sub match {
    my ($self, $seq) = @_;  ### we need to PASS IN THE PROFILE DB NAME SOMEHOW!!
    my @matches = ( );

    # Make sure that we received a sequence object.
    $self->throw("Not a Bio::Seq object")
	if((! ref($seq)) || (! $seq->isa("Bio::SeqI")));

    my $pf = $self->{profile}->{pfscan};
    my $profiles = $self->{profile}->{profiledb};
    my %localization = %{$self->{profile}->{profileids}};

    # Make a temp file for the sequence
    my $io = new Bio::Root::IO;
    my($fh, $tempfile) = $io->tempfile();
    my $out = Bio::SeqIO->new(-fh => $fh, -format => 'fasta');
    $out->write_seq($seq);
    $out->close();

    close($fh);

    # Open pfscan and use it as a filter
    open(RESULTS, "$pf -f $tempfile $profiles |") || $self->throw("Can't load pfscan: $!");
    

    # Read the results of pfscan
    my @seq = <RESULTS>;
    close RESULTS;

    #must parse the output
    #score    number pos.  x  - y  PROSITE_ID|DESCRIPTION

    # Is there any output from pfscan?
    if ($#seq >= 0) {
	
	# break each line down and extract what we need
	foreach my $line (@seq) { 
	    if($line=~/^\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+\-\s+(\d+)\s+(\S+)\|/) {
		my  $score = $1;
		my  $num = $2;
		my  $pos = $3;
		my  $start = $4;
		my  $end = $5;
		my  $profile_id = $6;
		my  $desc = $'; #'
		chomp ($desc);
		
=pod
 my %localization = ( PS50253 => "CytoplasmicMembrane", PS50283 => "CytoplasmicMembrane", PS50830 => "Extracellular", PS50847 => "Cell Wall", PS50850 => "CytoplasmicMembrane", PS50862 => "Cytoplasmic");
=cut


                 my $m = new Bio::Tools::PSort::Profile::Match(-profileid      => $profile_id,
							       -localization => $localization{$profile_id},
							       -comment      => $desc,
							       -start        => $start,
							       -end          => $end);
		 push(@matches, $m);
	    }#end if (regex)    
	 }#end foreach


	}#end if
   
  return @matches;

}




=head1 NAME

Bio::Tools::PSort::Profile - Perl implementation of the Profile protein subcellular
localization method.

=head1 SYNOPSIS

  use Bio::Tools::PSort::Profile;

  # Load a previously trained profile from a file.
  $motif = new Bio::Tools::PSort::Profile(-database => 'profile-db.txt');

  # Attempt to match on a Bio::Seq object.
  @matches = $motif->match($seq);
  print($seq->display_id . ": matched " . $_->profile_id . "\n") for(@matches);


=head1 DESCRIPTION

Bio::Tools::Motif uses a selection of motifs that have been identified
to be typical of proteins resident at a specific subcellular localization.
The module accepts a Bio::Seq object and attempts to match it against
a database, returning one or more Bio::Tools::Motif::Match objects with
the prediction information if successful.

=head1 CONSTRUCTOR

   $motif = new Bio::Tools::Motif(-database => 'motif-db.txt');

The Motif constructor accepts the name of an existing database file.

=head1 METHODS

   @matches = $motif->match($seq);

The match method accepts a Bio::Seq object as an argument and returns
an array of Bio::Tools::Motif::Match objects that matched the given
sequence.

 
=head1 AUTHOR

Fiona Brinkman, Chris Walsh, Matthew Laird <psort-mail@sfu.ca>
Brinkman Lab, Simon Fraser University, Canada

=head1 SEE ALSO

Bio::Tools::Motif::Pattern, Bio::Tools::Motif::Match

=head1 ACKNOWLEGEMENTS

Thanks go out to Fiona Brinkman, Jennifer Gardy and the other members of the
Simon Fraser University Brinkman laboratory.

=cut

1;
