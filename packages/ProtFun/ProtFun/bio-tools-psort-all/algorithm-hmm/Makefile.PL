use ExtUtils::MakeMaker;
use lib '../lib';
use Bio::Tools::PSort::Install;

my $CC = 'g++';

@libs = qw/hmmer squid/;
%paths = ();
$libpath = '-lm -lpthread';

$libpath .= Bio::Tools::PSort::Install->makeLibPath(@libs);

WriteMakefile('NAME'          => 'Algorithm::HMM',
	      'VERSION_FROM'  => 'lib/Algorithm/HMM.pm',
	      'PREREQ_PM'     => {},
	      ($] >= 5.005 ?
	       (ABSTRACT_FROM => 'lib/Algorithm/HMM.pm',
		AUTHOR        => 'Cory Spencer <cspencer@sfu.ca>') : ()),
	      'LIBS'          => $libpath,
	      'CC'            => $CC,
	      'LD'            => '$(CC)',
	      'INC'           => '-Isrc/include',
	      'OBJECT'        => "HMM.o hmm-binding.o",
	      'XSOPT'         => '-C++',
	      );
