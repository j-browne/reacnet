#include <string.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include "nucleus.h"
#include "globals.h"

size_t nucIndex(char* nuc) {
	size_t i;
	for (i=0; i<numNuclei; ++i) {
		if (strcmp(nuc,nuclei[i].name) == 0) {
			return i;
		}
	}
	return -1;
}

bool get_nuclei(const char* filename) {
	size_t maxNumNuclei=0;
	FILE* infile=fopen(filename,"r");
	if (infile==NULL) {
		fprintf(stderr,"Error opening %s.\n",filename);
		return false;
	}
	// Read in reactions
	while (!feof(infile)) {
		if (numNuclei>=maxNumNuclei) {
			if (maxNumNuclei<1) {
				maxNumNuclei=1;
			} else {
				maxNumNuclei*=2;
			}
			nuclei=realloc(nuclei,maxNumNuclei*sizeof(nucleus_t));
			abun=realloc(abun,maxNumNuclei*sizeof(double));
			if (nuclei==NULL) {
				fprintf(stderr,"Error allocating space for 'nuclei'.\n");
				return false;
			}
			if (abun==NULL) {
				fprintf(stderr,"Error allocating space for 'abun'.\n");
				return false;
			}
		}

		// Process file
		if (fscanf(infile, "%s\t%zu\t%zu\t%lf", nuclei[numNuclei].name, &nuclei[numNuclei].Z, &nuclei[numNuclei].N, &abun[numNuclei]) == 4) {
			++numNuclei;
		}
	}

	nuclei=realloc(nuclei,numNuclei*sizeof(nucleus_t));
	abun=realloc(abun,maxNumNuclei*sizeof(double));
	if (nuclei==NULL) {
		fprintf(stderr,"Error allocating space for 'nuclei'.\n");
		return false;
	}
	if (abun==NULL) {
		fprintf(stderr,"Error allocating space for 'abun'.\n");
		return false;
	}
	fclose(infile);
	return true;
}
