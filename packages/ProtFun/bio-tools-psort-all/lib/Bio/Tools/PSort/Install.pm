package Bio::Tools::PSort::Install;

use ExtUtils::Liblist;

sub findLib {
    my ($self, $lib) = @_;

    my @results = ExtUtils::Liblist->ext($lib, 0, 1);

    return $results[3];
}

sub findLD {
    my ($self) = @_;

    @paths = ();

    open(LD, "</etc/ld.so.conf") or die "Error opening /etc/ld.so.conf, this should never happen: $!\n";
    while(<LD>) {
	chomp;
	push @paths, $_ if -d $_;
    }
    close LD;

    if(-d '/etc/ld.so.conf.d/') {
	opendir($dh, '/etc/ld.so.conf.d/') || die "/etc/ld.so.conf.d/ exists but is not readable, that's weird: $!\n";
	while(readdir $dh) {
	    chomp;
	    push @paths, $self->checkDir($_) if -f $_;
	}

	close $dh;
    }

    return @paths;
}

sub checkDir {
    my ($self, $file) = @_;

    my @paths;

    open(LD, "<$file") or die "Error opening $file, that's weird... $!\n";
    while(<LD>) {
	chomp;
	push @paths, $_ if -d $_;
    }
    close LD;

    return @paths;
}

sub checkLib {
    my ($self, $lib, $sixtyfour) = @_;

    $libpath = '';
    $passed = 0;

    @ldpaths = Bio::Tools::PSort::Install->findLD();

    until($passed) {
	if($path = Bio::Tools::PSort::Install->findLib("-l" . $lib)) {
	    if(grep(/^$path$/, @ldpaths)) {
		$passed = 1;
	    } else {
		$libpath = $path;
		print "\nWe found $lib in $path, that's the good news.\n";
		print "The bad news is $lib is likely not in your dynamic linker path, you must do one of two things:\n";
		print "- have your system administrator add $path to '/etc/ld.so.conf' (RECOMMENDED!)\n";
		print "- add $path to the 'LD_LIBRARY_PATH' environment variable before each execution\n";
		print "\nPlease consult your system administrator or email us for further assistance.\n";
		$input = ExtUtils::MakeMaker::prompt("\nDo you want to continue anyways? [Y/n]", "Y");
		if($input =~ m/^Y$/i) {
		    $passed = 1;
		} else {
		    die "\nAlright exiting, you can try to run the install process again when you've corrected this problem\n";
		}
	    }
	} else {
	    $input = ExtUtils::MakeMaker::prompt("\n$lib was not found in the dynamic linker path, is there somewhere else we should check? ", "/usr/local/lib$sixtyfour");
	    if($path = Bio::Tools::PSort::Install->findLib("-L$input -l$lib")) {
		$passed = 1;
		$libpath = $input;
		
		if(grep(/^$path$/, @ldpaths)) {
 		    print "\nWe found $lib and $input is in your dynamic linker path, that's good.\n";
		    print "It's slightly odd that perl didn't find the library automatically\n";
		    print "however all should be fine.  If things break during 'make test' check your\n";
		    print "/etc/ld.so.conf to ensure $input is in it.\n";
  	        } else {
		    print "\nWe found $lib in $input, that's the good news.\n";
		    print "The bad news is $lib is likely not in your dynamic linker path, you must do one of two things:\n";
		    print "- have your system administrator add $input to '/etc/ld.so.conf' (RECOMMENDED!)\n";
		    print "- add $input to the 'LD_LIBRARY_PATH' environment variable before each execution\n";
		    print "\nPlease consult your system administrator or email us for further assistance.\n";
	        }
	    } else {
		$libpath = $input;
		print "\nWe're still not finding $lib either in your dynamic linker path (/etc/ld.so.conf) or\n";
		print "in $input.\n";
	    }
	    
	    $input = ExtUtils::MakeMaker::prompt("\nDo you want to continue anyways? [Y/n]", "Y");
	    if($input =~ m/^Y$/i) {
		$passed = 1;
	    } else {
		die "\nAlright exiting, you can try to run the install process again when you've corrected this problem\n";
	    }
	    
	}
    }

    return $libpath;
}

sub makeLibPath {
    ($self, @libs) = @_;

    $sixtyfour = '';
    my $libpath = '';
    %paths = ();

    if($^O =~ m/^linux$/i) {
	$sixtyfour = '64' if(-d '/usr/lib64');

	foreach my $lib (@libs) {
	    if($path = Bio::Tools::PSort::Install->checkLib($lib, $sixtyfour)) {
		$paths{$path} = 1;
	    }

	}

	foreach my $path (keys %paths) {
	    $libpath .= " -L$path";
	}

	foreach my $lib (@libs) {
	    $libpath .= " -l$lib";
	}
    } else {
	foreach my $lib (@libs) {
	    $libpath .= " -l$lib";
	}
    }

    return $libpath;
}

1;
