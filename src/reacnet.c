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

		double J[numNuclei][numNuclei];
		double b[numNuclei];
		double A[numNuclei*numNuclei];

		// Calculate Jacobian
		size_t j,k,l;
		for (j=0; j<numNuclei; ++j) {
			for (k=0; k<numNuclei; ++k) {
				J[j][j]=0;
			}
		}
		for (j=0; j<numReactions; ++j) {
			double R=reaction_rate(reactions[j],ops.temp);
			for (k=0; k<reactions[j].numIn; ++k) {
				R *= pow(ops.rho,reactions[j].numIn-1)*abun[reactions[j].in[k]];
			}
			for (k=0; k<reactions[j].numIn; ++k) {
				for (l=0; l<reactions[j].numIn; ++l) {
					J[reactions[j].in[l]][reactions[j].in[k]]=-R/abun[reactions[j].in[k]];
				}
				for (l=0; l<reactions[j].numOut; ++l) {
					J[reactions[j].out[l]][reactions[j].in[k]]=R/abun[reactions[j].in[k]];
				}
			}
		}

		// Solve (1/h-J)d=b for d
		// Calculate rhs (b)
		for (j=0; j<numNuclei; ++j) {
			b[j]=0;
		}
		for (j=0; j<numReactions; ++j) {
			double R=reaction_rate(reactions[j],ops.temp);
			for (k=0; k<reactions[j].numIn; ++k) {
				R *= pow(ops.rho,reactions[j].numIn-1)*abun[reactions[j].in[k]];
			}
			for (k=0; k<reactions[j].numIn; ++k) {
				b[reactions[j].in[k]]+=-R;
			}
			for (k=0; k<reactions[j].numOut; ++k) {
				b[reactions[j].out[k]]+=R;
			}
		}

		// Calculate lhs (1/h-J)
		// Also, reorder the array for FORTRAN
		for (j=0; j<numNuclei; ++j) {
			for (k=0; k<numNuclei; ++k) {
				A[k+numNuclei*j]=-J[k][j];
			}
			A[j+numNuclei*j] = A[j+numNuclei*j] + 1/ops.dt;
		}

		// Solve
		long c1=numNuclei;
		long c2=1;
		long pivot[numNuclei];
		long ok;
		dgesv_(&c1,&c2,A,&c1,pivot,b,&c1,&ok);

		// Update abundance
		for (j=0; j<numNuclei; ++j) {
			abun[j] += b[j];
		}

		// Output
		fprintf(outfile,"%le",t);
		print_abun(outfile);
		fprintf(outfile,"\n");
	}
	fclose(outfile);
}
