#ifndef OPTIONS_H
#define OPTIONS_H

#include <stdbool.h>
#include <stdint.h>

typedef struct
{
	char nuclei_filename[100];
	char reaclib_filename[100];
	char output_filename[100];
	uint32_t nSteps; // Number of steps in the simulation
	double dt; // Time step in s
	double temp; // Temperature in K
	double rho; // Density in g/cm^3
	bool seedFlag; // Whether to use a specific seed
	unsigned int seed;
} options_t;

bool get_ops(const char* filename);
#endif // OPTIONS_H
