#ifndef GLOBALS_H
#define GLOBALS_H

#include <stdlib.h>
#include "nucleus.h"
#include "reaction.h"
#include "options.h"

extern size_t numNuclei;
extern size_t numReactions;
extern nucleus_t* nuclei;
extern double* abun;
extern reaction_t* reactions;
extern options_t ops;

#endif // GLOBALS_H
