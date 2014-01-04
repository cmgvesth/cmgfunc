use ExtUtils::MakeMaker;
use lib '../lib';
use Bio::Tools::PSort::Install;

$CC = 'cc';

%args = ();

@libs = qw/modhmm/;
%paths = ();
$libpath = '-lm';

$libpath .= Bio::Tools::PSort::Install->makeLibPath(@libs);

WriteMakefile('NAME'          => 'Bio::Tools::PSort::ModHMM',
	      'VERSION_FROM'  => 'lib/Bio/Tools/PSort/ModHMM.pm',
	      'PREREQ_PM'     => {},
	      ($] >= 5.005 ?
	       (ABSTRACT_FROM => 'lib/Bio/Tools/PSort/ModHMM.pm',
		AUTHOR        => 'Matthew Laird <matt@brinkman.mbb.sfu.ca>') : ()),
	      'OPTIMIZE'      => '-O3',
	      'LIBS'          => $libpath,
	      'CC'            => $CC,
	      'LD'            => '$(CC)',
	      'OBJECT'        => 'modhmm.o binding.o',
	      'XSOPT'         => '',
	      %args);
