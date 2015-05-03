#include <stdbool.h>
#include <tgmath.h>
#include <stdio.h>
#include <string.h>
#include "reaction.h"
#include "nucleus.h"
#include "globals.h"

double reactionRate(reaction_t r, double T) {
	double exponent=r.param[0];
	size_t i;
	for (i=1; i<6; ++i) {
		exponent += r.param[i]*pow(T*1E-9, (2*i-5)/3.);
	}
	exponent += r.param[6]*log(T*1E-9);
	return exp(exponent);
}

bool get_reactions(const char* filename) {
	const int limits[11][2]={{1,1},{1,2},{1,3},{2,1},{2,2},{2,3},{2,4},{3,1},{3,2},{4,2},{1,4}}; // Structure of different reaclib chapters
	size_t maxNumReactions=0;
	FILE* infile=fopen(filename,"r");
	if (infile==NULL) {
		fprintf(stderr,"Error opening %s.\n",filename);
		return false;
	}
	// Read in reactions
	while (!feof(infile)) {
		int chap;
		char line[81];
		bool add=true;
		size_t i;

		if (numReactions>=maxNumReactions) {
			if (maxNumReactions<1) {
				maxNumReactions=1;
			} else {
				maxNumReactions*=2;
			}
			reactions=realloc(reactions,maxNumReactions*sizeof(reaction_t));
			if (reactions==NULL) {
				fprintf(stderr,"Error allocating space for 'reactions'.\n");
				return false;
			}
		}

		fscanf(infile,"%d",&chap);
		getc(infile);
		fgets(line,81,infile);

		if (chap<1 || chap>11) {
			add=false;
		} else {
			reactions[numReactions].numIn=limits[chap-1][0];
			reactions[numReactions].numOut=limits[chap-1][1];

			for (i=0;i<(reactions[numReactions].numIn+reactions[numReactions].numOut);++i) {
				char name[6];
				char* namestart=name;
				strncpy(name,line+5*(i+1),5);
				name[5]='\0';
				while (*namestart == ' ') {
					++namestart;
				}
				if (i<reactions[numReactions].numIn) {
					if ((reactions[numReactions].in[i]=nucIndex(namestart,nuclei,numNuclei))==-1) {
						add=false;
					}
				} else {
					if ((reactions[numReactions].out[i-reactions[numReactions].numIn]=nucIndex(namestart,nuclei,numNuclei))==-1) {
						add=false;
					}
				}
			}

			strncpy(reactions[numReactions].label,line+43,4); // TODO: Get rid of spaces
			reactions[numReactions].label[4]='\0';
			reactions[numReactions].res=line[47];
			reactions[numReactions].rev=(line[48]=='v');
			sscanf(line+52,"%12le",&reactions[numReactions].Q);
		}
		for (i=0;i<7;++i) {
			fscanf(infile,"%13le",&reactions[numReactions].param[i]);
		}

		if (add) {
			++numReactions;
		}
	}

	reactions=realloc(reactions,numReactions*sizeof(reaction_t));
	if (reactions==NULL) {
		fprintf(stderr,"Error allocating space for 'reactions'.\n");
		return false;
	}
	fclose(infile);
	return true;
}

