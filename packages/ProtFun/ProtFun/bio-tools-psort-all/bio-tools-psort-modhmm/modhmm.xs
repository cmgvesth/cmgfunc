#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#include "modhmm0.92b/structs.h"

MODULE = Bio::Tools::PSort::ModHMM	PACKAGE = Bio::Tools::PSort::ModHMM

int
find_helices(seq, hmmfile, repfile, path)
	char *seq
	char *hmmfile
	char *repfile
	char *path
   OUTPUT:
	RETVAL

Bio::Tools::PSort::ModHMM::HelixSites *
get_helices(seq, hmmfile, repfile, path)
	char *seq
	char *hmmfile
	char *repfile
	char *path
   OUTPUT:
	RETVAL


MODULE = Bio::Tools::PSort::ModHMM	PACKAGE = Bio::Tools::PSort::ModHMM::HelixSitesPtr	PREFIX = helix_

PROTOTYPES: ENABLE

int
helix_count(hSites)
	Bio::Tools::PSort::ModHMM::HelixSites *	hSites

int
helix_start(hSites, helix)
	Bio::Tools::PSort::ModHMM::HelixSites *	hSites
	int helix

int
helix_end(hSites, helix)
	Bio::Tools::PSort::ModHMM::HelixSites *	hSites
	int helix

void
helix_DESTROY(hSites)
	Bio::Tools::PSort::ModHMM::HelixSites * hSites
