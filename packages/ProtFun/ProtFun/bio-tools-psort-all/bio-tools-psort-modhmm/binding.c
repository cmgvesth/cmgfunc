#include <stdio.h>

#include "modhmm0.92b/structs.h"

int helix_count(helix_sites * hSites)
{
    if(hSites) {
	return hSites->helix_count;
    }

    return -1;
}

int helix_start(helix_sites * hSites, int helix)
{
    if(!hSites) {
	return -1;
    }

    // We're doing 1 based indexing for the class users
    // but 0 based internally
    helix--;

    if(helix <= hSites->helix_count) {
	return hSites->helix[helix].start;
    }

    return -1;
}

int helix_end(helix_sites * hSites, int helix)
{
    if(!hSites) {
	return -1;
    }

    // We're doing 1 based indexing for the class users
    // but 0 based internally
    helix--;

    if(helix <= hSites->helix_count) {
	return hSites->helix[helix].end;
    }

    return -1;
}

void helix_DESTROY(helix_sites * hSites)
{
    if(hSites) {
	if(hSites->helix) {
	    free(hSites->helix);
	}

	free(hSites);
    }
}
