#ifndef REACTION_H
#define REACTION_H

#include <stdbool.h>

typedef struct
{
	size_t numIn;
	size_t numOut;
	short in[4];
	short out[4];
	double param[7];
	char label[5];
	char res; // Non-resonant=' ' or 'n', Resonant='r', Weak='w'
	bool rev;
	double Q;
} reaction_t;

double reaction_rate(reaction_t r, double T);
bool get_reactions(const char* filename);

#endif // REACTION_H
