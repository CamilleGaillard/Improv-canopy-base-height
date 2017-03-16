#ifndef SeedBankClassConstants_h___
#define SeedBankClassConstants_h___
#include "PlantPopClassConstants.h"

const int    MAX_SEED_NUM       = 40;
int const    MAX_SEED_BANK_SIZE = MAX_POP_SIZE*MAX_SEED_NUM; //40*40*100;

#ifdef FIXED_TRAITS

double const MUT_PROBABILITY    = -10.0;
double const MUT_RATE           = 0.;
double const MUT_RATE_LO        = 1. - MUT_RATE;
double const MUT_RATE_UP        = 1. + MUT_RATE;
double const CROSSOVER_CR       = -10.;
double const CROSSOVER_F        = 0.;

#else  // FIXED_TRAITS not def

double const MUT_PROBABILITY    = 0.01;
double const MUT_RATE           = 0.05;
double const MUT_RATE_LO        = 1. - MUT_RATE;
double const MUT_RATE_UP        = 1. + MUT_RATE;

double const CROSSOVER_CR       = 0.5;  // crossover probability (for single traits)

// not used

double const CROSSOVER_F        = 0.5;  // crossover rate

double const CROSSOVER_ALL      = 0.2;
double const CROSSOVER_PROB     = 0.5;   // crossover in new model

double const CROSSOVER_EXT      = 0.05;  // crossover with other species if no appropriate seed found

// double const MUT_PROBABILITY    = 0.0;
// double const MUT_RATE           = 0.;
// double const MUT_RATE_LO        = 1. - MUT_RATE;
// double const MUT_RATE_UP        = 1. + MUT_RATE;
// double const CROSSOVER_CR       = 0.;
// double const CROSSOVER_F        = 0.;

#endif  // FIXED_TRAITS


#endif









