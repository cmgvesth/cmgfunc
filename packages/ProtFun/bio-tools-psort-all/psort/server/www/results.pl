#! /usr/bin/perl -w

use Bio::Tools::PSort::Report::Formatter;
use Bio::Tools::PSort::XMLRPC::Client;
use Bio::SeqIO;
use Bio::Seq;

use IO::String;
use CGI;

use strict;

my $RPCSERVER = "http://localhost/RPC";
my $MAXFSIZE  = 300000;
my $CUTOFF    = 0.75;


MAIN: {
  my ($cgi, @errors);
  my $output;
  my @seqs;

  $cgi = new CGI;

  # Check to see if something's been submitted to use.
  if($cgi->param('submit')) {
    # Get the parameters sent to us.
    my $fname = $cgi->param('filename');
    my $seqs  = $cgi->param('seqs');
    my $fmt   = $cgi->param('format') || 'html';
    my $gram  = $cgi->param('gram') || 'none';
    push(@errors, "You have not selected if the sequences are Gram positive or negative, please press the back button and select the correct Gram.") if ($gram eq 'none');
    my $analysis = ($gram eq 'positive'? 'psortp' : 'psortn');

    my $output = ($gram eq 'positive'? 'bayesp' : 'bayesn');
    my $mcutoff = ($gram eq 'positive'? 0.5 : 0.4);
    my $ogram = ($gram eq 'positive'? ' grampos' : 'gramneg');
    my @amods = ($gram eq 'positive'? 
		 qw(CMSVM+ CWSVM+ CytoSVM+ ECSVM+ HMMTOP+ Motif+ Profile+ SCL-BLAST+ SCL-BLASTe+ Signal+) :
		 qw(CMSVM- CytoSVM- ECSVM- HMMTOP- Motif- OMPMotif- OMSVM- PPSVM- SCL-BLAST- SCL-BLASTe- Signal-));

    my $fh = $cgi->upload('filename');
    if($fname) {
      if((-s $fname) < $MAXFSIZE) {
	$seqs = join('', <$fh>);
      } else {
	push(@errors, "The file you have uploaded exceeds the maximum allowable size of $MAXFSIZE bytes.  Please split your sequences up into smaller files, and upload them individually.");
      }
    }

    push(@errors, "The number of sequences you have submitted exceeds the maximum allowable size of $MAXFSIZE bytes.  Please split your sequences into smaller sets, and upload them individually.")
    if(length($seqs) > $MAXFSIZE);

    if(! @errors) {
      # If no sequence ID's were provided in the form, then automagically assign
      # our own.
      my $id = 1;
      # Convert mac format to unix
      unless ($seqs =~ /\n/) {
          $seqs =~ s/\r/\n/g;
      }

      $seqs =~ s/>\s*\n/sprintf(">Unknown_%d\n", $id++)/ge;

      # Read in all the sequences.
      eval {
	my $seqfh = new IO::String($seqs);
	my $seqio = new Bio::SeqIO('-fh' => $seqfh, '-format' => 'fasta');
	my $invalid;
	while(my $seq = $seqio->next_seq) {
	  $invalid++ if(! $seq->validate_seq());
	  push(@seqs, $seq);
	}

	# If of counter indicated that validation for one or more sequences
	# failed, then set an error message.
	push(@errors, "One or more sequences were in an invalid format.  Please ensure that your sequences are in FASTA format and that the sequence data contains only valid characters.")
	if($invalid);
      };
      push(@errors, "An error occured while processing the sequence stream.") if($@);

      # If we've got some sequences, then let's PSort!
      if(@seqs && (! @errors)) {
	my ($psort, @res, $fmtr);

	# Create a new PSort object.
	$psort = new Bio::Tools::PSort::XMLRPC::Client(-server => $RPCSERVER);

	# Push all the reports onto an array.
	push(@res, $psort->classify($_, { analysis => $analysis, output => $output })) for(@seqs);

	# Create a new report formatter of the specified type.
	$fmtr   = new Bio::Tools::PSort::Report::Formatter();
	$output = $fmtr->format($fmt, {-gram => $ogram, -amods => \@amods, -mcutoff => $mcutoff }, @res);
      } elsif(! @seqs) {
	push(@errors, "No sequences were submitted.");
      }
    }

    print($cgi->header('text/plain'));
    if(@errors) {
      $output = "One or more errors occurred during your request:<br>"
	        . "<ul><li>" . join('<li>', @errors) . "</ul>";
    }
    print($output);
  }
}
