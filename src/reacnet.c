#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <tgmath.h>
#include <time.h>
#include "nucleus.h"
#include "reaction.h"
#include "options.h"
#include "globals.h"
#include "lapack.h"

size_t numNuclei=0;
size_t numReactions=0;
nucleus_t* nuclei=NULL;
double* abun=NULL;
reaction_t* reactions=NULL;
options_t ops;

void print_abun(FILE* stream);
void sim();

int main (int argc, char* argv[]) {
	// Initialize options
	int i;
	get_ops(NULL);
	for (i=1; i<argc; ++i) {
		get_ops(argv[i]);
	}

	// Seed RNG
	if (ops.seedFlag) {
		srand(ops.seed);
	} else {
		srand(time(NULL));
	}

	// Read in nuclei
	if (get_nuclei(ops.nuclei_filename)==false) {
		exit(EXIT_FAILURE);
	}

	// Read in reactions
	if (get_reactions(ops.reaclib_filename)==false) {
		exit(EXIT_FAILURE);
	}

	sim();

	free(nuclei);
	free(abun);
	free(reactions);
	return 0;
}

void print_abun(FILE* stream) {
	size_t j;
	for (j=0; j<numNuclei; ++j) {
		fprintf(stream,"\t%e",abun[j]);
	}
	return;
}

void sim() {
	// Header for output
	FILE* outfile=fopen(ops.output_filename,"w");
	if (outfile==NULL) {
		fprintf(stderr,"Error opening %s.\n",ops.output_filename);
		exit(EXIT_FAILURE);
	}
	size_t i;
	fprintf(outfile,"# time");
	for (i=0; i<numNuclei; ++i) {
		fprintf(outfile,"\t%s",nuclei[i].name);
	}
	fprintf(outfile,"\n");

	double t=0;
	for (i=0; i<ops.nSteps; ++i) {
		t += ops.dt;
		fprintf(outfile,"%lf",t);
		print_abun(outfile);
		fprintf(outfile,"\n");
	}
	fclose(outfile);
}
