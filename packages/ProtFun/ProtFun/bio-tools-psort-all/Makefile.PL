use ExtUtils::MakeMaker;
use lib 'lib';
use Bio::Tools::PSort::Install;
use File::Spec;
use Cwd;
use File::Copy;
#use strict;

my %packages = (
	     'XML::RPC'     => '0.0/XML Based RPC client/making calls to remote PSortb servers, if you don\'t plan to do this, ignore this warning/XML::RPC',
	     );

my $BLASTDIR = '';
my $PFTOOLS = '';
my $PSORT_ROOT = '';

# Borrowed from Bioperl
sub check_package {
    my ($name,$str) = @_;
    my ($ver,$desc,$expl,$module) = split( /\// , $str);
    if( !eval "require($name)") {
	print "External Module $name, $desc, is not installed on this computer.\n  The $module in PSortb needs it for $expl\n\n";
	return 1;
    }

    return 0;
}

sub get_path_info {
    if(!eval "require Bio::Root::IO") {
	die "Error, Bioperl doesn't seem to be installed!\n";
    }

    # First let's check for BLAST
    find_blast();

    # Do we have pfscan?
    $PFTOOLS = find_file('pfscan', 0);

    get_psort_root();
}

sub find_blast {

    if(exists($ENV{BLASTDIR})) {
	# Good, the user set it up for us!
	$BLASTDIR = $ENV{BLASTDIR};
	$BLASTDIR = prompt("We think we've found blastall, please enter the path if this isn't correct", "$BLASTDIR");
    } elsif(my $str = Bio::Root::IO->exists_exe('blastall')) {
	my ($volume,$directories,$file) = File::Spec->splitpath($str);
	$BLASTDIR = prompt("We think we've found blastall, please enter the path if this isn't correct", "$directories");
    } else {
	$BLASTDIR = prompt("We couldn't find blastall, please enter the path", '');
    }

    return 1
	if( (-d $BLASTDIR) && ( -f "$BLASTDIR/blastall" ) );

    if( (-d "$BLASTDIR/bin") && ( -f "$BLASTDIR/bin/blastall" ) ) {
	$BLASTDIR = "$BLASTDIR/bin";
	my $ask = prompt("We think we've found blastall in $BLASTDIR, is that were the executable is?", 'Y');
	if($ask =~ /[Yy]/) {
	    return 1;
	}
    }

    die "I'm sorry, $BLASTDIR doesn't seem to contain blastall, we don't know where blastall is, aborting\n";

}

sub find_file {
    my ($exe, $required) = @_;
    my $newpath = '';

    if(my $str = Bio::Root::IO->exists_exe($exe)) {
	my ($volume,$directories,$file) = File::Spec->splitpath($str);	
	$newpath = prompt("We think we've found $exe, please enter the path if this isn't correct", "$directories");
    } else {
	if($required) {
	    $newpath = prompt("We couldn't find $exe, please enter the path", '');
	} else {
	    print "We couldn't find $exe, is it installed?  While this program isn't strictly required, ";
	    print "it's strongly recommended you stop now and install it first.  Missing modules will ";
	    print "slightly compromise accuracy of PSortb.";
	    $newpath = prompt(" Please enter the path for $exe.  If you don't have $exe installed and wish to disable this module please simply hit enter.", '');	    
	}
    }
    
    return $newpath
	if( (-d $newpath) && ( -f "$newpath/$exe" ) );

    if($required) {
	die "I'm sorry, $newpath doesn't seem to contain $exe, we don't know where $exe is, aborting\n";
    } else {
	return "NOTINSTALLED";
    }
}

sub get_psort_root {
    my $root = prompt("Where do you plan to install the psort configuration files?", '/usr/local/psortb');
    $PSORT_ROOT = $root;
}

sub update_psortb {
    print "Updating psort/bin/psort with path hints\n";

    # Update psortb script here
    unless (open(PSORT, "<psort/bin/psort.tmpt")) {
	print "Error: We weren't able to open the psort script which should be in psort/bin/psort: $!\n";
	print "We can't set the path hints in the script, you will have to set the PSORT_ROOT, PFTOOLS_ROOT, and BLASTDIR environmental variables\n";
	return;
    }

    my @psort = ();
    while(<PSORT>) {
#	chomp;
#	print $_;
	push @psort, $_;
    }

#    my @psort = <PSORT>;

    close PSORT;

    unless(open(PSORT, ">psort/bin/psort")) {
	print "Error: We weren't able to open the psort script located in psort/bin/psort for writing: $!\n";
	print "We can't set the path hints in the script, you will have to set the PSORT_ROOT, PFTOOLS_ROOT, and BLASTDIR environmental variables\n";
	return;
    }

    foreach my $line (@psort) {
	if($line =~ /^\s*my \$root = '';/) {
	    print PSORT "  my \$root = '$PSORT_ROOT';\n";
	} elsif($line =~ /^\s*my \$blastdir = '';/) {
	    print PSORT "  my \$blastdir = '$BLASTDIR';\n";
	} elsif($line =~ /^\s*my \$pftools = '';/) {
	    print PSORT "  my \$pftools = '$PFTOOLS';\n";
	} else {
	    print PSORT $line;
	}
    }

    close PSORT;
}

sub check_psort_path {
    my $pwd = &Cwd::cwd();
    if(($pwd =~ /^$PSORT_ROOT/) || ($PSORT_ROOT =~ /^$pwd/)) {
	print "$PSORT_ROOT matches our current directory, do you plan to run it out of the current directory?\n";
	print "That's ok, we're just pointing it out.\n";
	return;
    } elsif( -e $PSORT_ROOT ) {
	print "Warning, $PSORT_ROOT exists, are you reinstalling or upgrading Psortb?\n";
	print "That's ok, we're just pointing it out.\n";
	return;
    }

#    print "Ok, skipping install, you will still have to copy the psort directory to where you want it to live.\n";
}

sub make_sclblast {

    print "Making sclblast databases\n";

    unless( -f "$BLASTDIR/formatdb" ){
	print "Error: We weren't able to find formatdb (part of BLAST) in $BLASTDIR,\n";
	print "you will need to manually run makedb.sh under psort/conf/analysis/sclblast/\n";
	return;
    }

    unless( -d 'psort/conf/analysis/sclblast/' ) {
	print "Error: We can't seem to find psort/conf/analysis/sclblast/ in the current directory,\n";
	print "perhaps you're not in the bio-tools-psort-all directory?\n";
	print "You will need to manually run makedb.sh under psort/conf/analysis/sclblast/\n";
	return;
    }

    chdir "psort/conf/analysis/sclblast/gramneg";
    system("$BLASTDIR/formatdb -i sclblast -p T -o T");
    chdir "../grampos";
    system("$BLASTDIR/formatdb -i sclblast -p T -o T");
    chdir "../archaea";
    system("$BLASTDIR/formatdb -i sclblast -p T -o T");
    chdir "../../../../../";
}

# Actual executing code...

foreach my $name (keys %packages) {
    check_package($name, $packages{$name});
}

get_path_info();

print "\n";

WriteMakefile('NAME'          => 'Bio::Tools::PSort',
	      'VERSION_FROM'  => 'lib/Bio/Tools/PSort.pm',
	      'PREREQ_PM'     => {'IO::String'     => 0},
	      ($] >= 5.005 ?
	       (ABSTRACT_FROM => 'lib/Bio/Tools/PSort.pm',
		AUTHOR        => 'Matthew Laird <lairdm@sfu.ca>') : ()),
	      DIR             => ['algorithm-hmm',
				  'bio-tools-psort-svmloc',
				  'bio-tools-psort-modhmm'],
	      );

check_psort_path();
update_psortb();
make_sclblast();

print "\nReady for make, next steps you must run:\n";
print "\$ make\n";
print "\$ make test\n";
print "\$ make install\n";
print "RECOMMENDED: \n";
print "cp -r psort $PSORT_ROOT\n";
print "\nPSortb will then be available at $PSORT_ROOT/bin/psort or in your path depending on your configuration\n\n";
