#ifndef NUCLEUS_H
#define NUCLEUS_H

#include <stddef.h>
#include <stdbool.h>

typedef struct
{
	char name[6];
	size_t Z;
	size_t N;
} nucleus_t;

size_t nucIndex(char* nuc, nucleus_t* nuclei, size_t numNuclei);
bool get_nuclei(const char* filename);

#endif // NUCLEUS_H
