#include <stdbool.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include "options.h"
#include "globals.h"

bool get_ops(const char* filename) {
	// Defaults
	if (filename==NULL) {
		strncpy(ops.nuclei_filename,"nuclei",sizeof(ops.nuclei_filename));
		strncpy(ops.reaclib_filename,"reaclib",sizeof(ops.reaclib_filename));
		strncpy(ops.output_filename,"out",sizeof(ops.output_filename));
		ops.nSteps=0;
		ops.dt=1e-15;
		ops.temp=1e8;
		ops.rho=1.622e2;
		ops.seedFlag=false;

		return true;
	}

	FILE* infile=fopen(filename,"r");
	if (infile==NULL) {
		fprintf(stderr,"Error opening %s.\n",filename);
		return false;
	}
	// Read in options
	while (!feof(infile)) {
		char key[100];
		char val[100];

		fscanf(infile,"%s\t%s",key,val);
		if (strcmp(key,"nuclei_filename") == 0) {
			strncpy(ops.nuclei_filename,key,sizeof(ops.nuclei_filename));
		} else if (strcmp(key,"reaclib_filename") == 0) {
			strncpy(ops.reaclib_filename,key,sizeof(ops.reaclib_filename));
		} else if (strcmp(key,"output_filename") == 0) {
			strncpy(ops.output_filename,key,sizeof(ops.output_filename));
		} else if (strcmp(key,"nSteps") == 0) {
			ops.dt=atoi(key);
		} else if (strcmp(key,"dt") == 0) {
			ops.dt=atof(key);
		} else if (strcmp(key,"temp") == 0) {
			ops.temp=atof(key);
		} else if (strcmp(key,"rho") == 0) {
			ops.rho=atof(key);
		} else if (strcmp(key,"seed") == 0) {
			ops.seedFlag=true;
			ops.seed=atoi(key);
		}
	}
	return true;
}
