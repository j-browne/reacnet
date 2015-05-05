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
	double t=0;
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
	fprintf(outfile,"%le",t);
	print_abun(outfile);
	fprintf(outfile,"\n");

	for (i=0; i<ops.nSteps; ++i) {
		t+=ops.dt;

		size_t j,k;
		double b[numNuclei];
		for (j=0; j<numNuclei; ++j) {
			b[j]=0;
		}
		for (j=0; j<numReactions; ++j) {
			// Determine rate
			double R=pow(ops.rho,reactions[j].numIn-1)*reaction_rate(reactions[j], ops.temp);
			for (k=0; k<reactions[j].numIn; ++k) {
				R*=abun[reactions[j].in[k]]/(nuclei[reactions[j].in[k]].Z+nuclei[reactions[j].in[k]].N);
			}
			// Use up reactants
			for (k=0; k<reactions[j].numIn; ++k) {
				b[reactions[j].in[k]]-=R*ops.dt*(nuclei[reactions[j].in[k]].Z+nuclei[reactions[j].in[k]].N);
			}
			// Create products
			for (k=0; k<reactions[j].numOut; ++k) {
				b[reactions[j].out[k]]+=R*ops.dt*(nuclei[reactions[j].in[k]].Z+nuclei[reactions[j].in[k]].N);
			}
		}

		for (j=0; j<numNuclei; ++j) {
			abun[j]+=b[j];
			if (abun[j]<0) {
				abun[j]=0;
			}
		}

		fprintf(outfile,"%e",t);
		print_abun(outfile);
		fprintf(outfile,"\n");
	}
	fclose(outfile);
}
