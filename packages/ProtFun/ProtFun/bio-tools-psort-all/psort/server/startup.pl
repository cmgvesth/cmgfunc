#! /usr/bin/perl

use Apache::Bio::Tools::PSort ();
use Apache::RPC::Server ();
use Apache::Registry ();
use Apache::Constants ();
use Apache ();

use strict;

my $PSORT_ROOT = "/usr/local/psort/";
my $BLAST_DIR  = "/path/to/BLAST/installation/";
my $PFTOOLS_DIR = "/path/to/PFTOOLS/installation/";

our %ENV;

MAIN: {
  my $hmmfile = "$PSORT_ROOT/conf/analysis/modhmm/S_TMHMM_0.92b.hmg";
  my $repfile = "$PSORT_ROOT/conf/analysis/modhmm/replacement_letter_multi.rpl";
  my $hmmpath = "$PSORT_ROOT/conf/analysis/modhmm/";
  my $motifdbp = "$PSORT_ROOT/conf/analysis/motif/grampos/motifs.txt";
  my $motifdbn = "$PSORT_ROOT/conf/analysis/motif/gramneg/motifs.txt";
  my $sighmmp  = "$PSORT_ROOT/conf/analysis/signal/grampos/model.hmm";
  my $sigsvmp  = "$PSORT_ROOT/conf/analysis/signal/grampos/model.svm";
  my $sighmmn  = "$PSORT_ROOT/conf/analysis/signal/gramneg/model.hmm";
  my $sigsvmn  = "$PSORT_ROOT/conf/analysis/signal/gramneg/model.svm";
  my $sigprogp = "$PSORT_ROOT/conf/analysis/signal/grampos/check-sig";
  my $sigprogn = "$PSORT_ROOT/conf/analysis/signal/gramneg/check-sig";
  my $scldbp   = "$PSORT_ROOT/conf/analysis/sclblast/grampos/sclblast";
  my $scldbn   = "$PSORT_ROOT/conf/analysis/sclblast/gramneg/sclblast";
  my $ompdb    = "$PSORT_ROOT/conf/analysis/omp-motif/omp-motifs.txt";
  my $pfscan   = "$PFTOOLS_DIR/pfscan";
  my $profiledbp  = "$PSORT_ROOT/conf/analysis/profile/grampos/ps_ALL";
  my $profileidsp = "$PSORT_ROOT/conf/analysis/profile/grampos/profile_ids";
  my $profiledbn  = "$PSORT_ROOT/conf/analysis/profile/gramneg/ps_ALL";
  my $profileidsn = "$PSORT_ROOT/conf/analysis/profile/gramneg/profile_ids";

  my $slmodelCn    = "$PSORT_ROOT/conf/analysis/subloc/gramneg/Cytoplasmic/SVM_MODEL.txt";
  my $slpatternsCn = "$PSORT_ROOT/conf/analysis/subloc/gramneg/Cytoplasmic/fre_patterns.txt";
  my $slmodelIn    = "$PSORT_ROOT/conf/analysis/subloc/gramneg/Innermembrane/SVM_MODEL.txt";
  my $slpatternsIn = "$PSORT_ROOT/conf/analysis/subloc/gramneg/Innermembrane/fre_patterns.txt";
  my $slmodelOn    = "$PSORT_ROOT/conf/analysis/subloc/gramneg/Outermembrane/SVM_MODEL.txt";
  my $slpatternsOn = "$PSORT_ROOT/conf/analysis/subloc/gramneg/Outermembrane/fre_patterns.txt";
  my $slmodelEn    = "$PSORT_ROOT/conf/analysis/subloc/gramneg/Extracellular/SVM_MODEL.txt";
  my $slpatternsEn = "$PSORT_ROOT/conf/analysis/subloc/gramneg/Extracellular/fre_patterns.txt";
  my $slmodelPn     = "$PSORT_ROOT/conf/analysis/subloc/gramneg/Periplasmic/SVM_MODEL.txt";
  my $slpatternsPn  = "$PSORT_ROOT/conf/analysis/subloc/gramneg/Periplasmic/fre_patterns.txt";

  my $slmodelCp     = "$PSORT_ROOT/conf/analysis/subloc/grampos/Cytoplasmic/SVM_MODEL.txt";
  my $slpatternsCp  = "$PSORT_ROOT/conf/analysis/subloc/grampos/Cytoplasmic/fre_patterns.txt";
  my $slmodelMp     = "$PSORT_ROOT/conf/analysis/subloc/grampos/Membrane/SVM_MODEL.txt";
  my $slpatternsMp  = "$PSORT_ROOT/conf/analysis/subloc/grampos/Membrane/fre_patterns.txt";
  my $slmodelWp     = "$PSORT_ROOT/conf/analysis/subloc/grampos/Cellwall/SVM_MODEL.txt";
  my $slpatternsWp  = "$PSORT_ROOT/conf/analysis/subloc/grampos/Cellwall/fre_patterns.txt";
  my $slmodelEp     = "$PSORT_ROOT/conf/analysis/subloc/grampos/Extracellular/SVM_MODEL.txt";
  my $slpatternsEp  = "$PSORT_ROOT/conf/analysis/subloc/grampos/Extracellular/fre_patterns.txt";

  my $bayesp   = "$PSORT_ROOT/conf/output/bayesian/grampos/bayes.model";
  my $bayesn   = "$PSORT_ROOT/conf/output/bayesian/gramneg/bayes.model";

  $ENV{BLASTDIR} = $BLAST_DIR;
  $ENV{PSORT_PFTOOLS} = $PFTOOLS_DIR;

  my $psort = Apache::Bio::Tools::PSort->instance();

  # Load the various modules.  (Note that the second parameter is the modules
  # "aliased" name (to allow more than one instance of a given module,
  # possibly with different params, in a pathway) and **MUST** agree with the
  # names in the Bayes model file.  **IF THEY DON'T AGREE YOUR RESULTS WILL
  # BE BUNG!**

  # Comparison function for SCL-BLASTe
  my $SCLSub = sub {
      my $res = shift;
      
      if($res->localization =~ /Unknown/) {
	  return 1;
      }
      
      return 0;
  };

  # Gram negative
  $psort->install('SCLBlast', 'SCL-BLAST-', { -database => $scldbn });
  $psort->install('SCLBlast', 'SCL-BLASTe-', { -database => $scldbn,
					       -exact => 1 });
  $psort->install('ModHMM', 'ModHMM-', { -hmmfile => $hmmfile,
					 -repfile => $repfile,
					 -path => $hmmpath,
					 -loc => 'CytoplasmicMembrane' });
  $psort->install('OMPMotif', 'OMPMotif-', { -database => $ompdb });
  $psort->install('Motif', 'Motif-', { -database => $motifdbn });
  $psort->install('Profile', 'Profile-', {-database => $profiledbn,
					 -program => $pfscan, 
					 -profileids => $profileidsn } );
  $psort->install('SVMLoc', 'CytoSVM-', { -model => $slmodelCn, 
					 -patterns => $slpatternsCn});
  $psort->install('SVMLoc', 'ECSVM-', { -model => $slmodelEn, 
					 -patterns => $slpatternsEn});
  $psort->install('SVMLoc', 'CMSVM-', { -model => $slmodelIn, 
					 -patterns => $slpatternsIn});
  $psort->install('SVMLoc', 'OMSVM-', { -model => $slmodelOn, 
					 -patterns => $slpatternsOn});
  $psort->install('SVMLoc', 'PPSVM-', { -model => $slmodelPn, 
					       -patterns => $slpatternsPn});
  $psort->install('Signal', 'Signal-', { -svm => $sigsvmn, -hmm => $sighmmn,
					-program => $sigprogn });
  $psort->add_path("psortn", "analysis", ['SCL-BLASTe-', $SCLSub,
					  ['SCL-BLAST-', sub { 1 },
					   ['ModHMM-', sub { 1 },
					    ['OMPMotif-', sub { 1 },
					     ['Motif-', sub { 1 },
					      ['Profile-', sub { 1 },
					       ['CytoSVM-', sub { 1 } ,
						['CMSVM-', sub { 1 } ,
						 ['PPSVM-', sub { 1 } ,
						  ['OMSVM-', sub { 1 } ,
						   ['ECSVM-', sub { 1 } ,
						    ['Signal-']]]]]]]]]]]]);

  $psort->install('Bayesian', 'Bayesian-', {-model => $bayesn });
  $psort->add_path("bayesn", "output", ['Bayesian-', sub { 1 }]);

  # Gram positive
  $psort->install('ModHMM', 'ModHMM+', { -hmmfile => $hmmfile,
					 -repfile => $repfile,
					 -path => $hmmpath,
					 -loc => 'CytoplasmicMembrane' });
  $psort->install('SCLBlast', 'SCL-BLAST+', { -database => $scldbp });
  $psort->install('SCLBlast', 'SCL-BLASTe+', { -database => $scldbp,
					       -exact => 1 });
  $psort->install('Motif', 'Motif+', { -database => $motifdbp });
  $psort->install('Profile', 'Profile+', {-database => $profiledbp,
					 -program => $pfscan, 
					 -profileids => $profileidsp } );
  $psort->install('SVMLoc', 'CytoSVM+', { -model => $slmodelCp, 
					 -patterns => $slpatternsCp});
  $psort->install('SVMLoc', 'ECSVM+', { -model => $slmodelEp, 
					 -patterns => $slpatternsEp});
  $psort->install('SVMLoc', 'CMSVM+', { -model => $slmodelMp, 
					 -patterns => $slpatternsMp});
  $psort->install('SVMLoc', 'CWSVM+', { -model => $slmodelWp, 
					 -patterns => $slpatternsWp});
  $psort->install('Signal', 'Signal+', { -svm => $sigsvmp, -hmm => $sighmmp,
					-program => $sigprogp, -gram => 'Positive' });
  $psort->add_path("psortp", "analysis", ['SCL-BLASTe+', $SCLSub,
					  ['ModHMM+', sub { 1 },
					   ['SCL-BLAST+', sub { 1 },
					    ['Motif+', sub { 1 },
					     ['Profile+', sub { 1 },
					      ['CytoSVM+', sub { 1 } ,
					       ['CMSVM+', sub {1} ,
						['ECSVM+', sub {1},
						 ['CWSVM+', sub {1},
						  ['Signal+']
						  ]]]]]]]]]);

  $psort->install('Bayesian', 'Bayesian+', {-model => $bayesp });
  $psort->add_path("bayesp", "output", ['Bayesian+', sub { 1 }]);

  return 1;
}
