#ifndef ParamFitDataStructures
#define ParamFitDataStructures

#include "Constants.h"

#include "ChargeFitDS.h"

struct ParameterFit {
  int nconf;                // The number of conformations in this fitting set
  int nqseed;               // The number of charge values to try for each
                            //   variable group
  int nqvar;                // The number of variable charge groups
  int nsigvar;              // The number of variable sigma groups
  int maxiter;              // The number of parameter optimization cycles
  int* qseed;               // The charge seeding counters (ticks forward with
                            //   each successive charge optimization trial)
  double qspc;              // The spacing of different charge trial starting
                            //   positions in all variable groups
  double wttemp;            // Temperature for assigning Boltzmann weights to
                            //   each conformation
  double* qval;             // The values of charges in this optimized set
  double* bqval;            // The best values of charges in any optimized set
                            //   thus far
  nail* qvar;               // Charge variability terms, one for each group
  nail* sigvar;             // Lennard-Jones sigma variability terms
  char sysconf[MAXNAME];    // File containing all system conformations
  char sysnrg[MAXNAME];     // File containing energies of each conformation
  char syswt[MAXNAME];      // Optional file containing weights of each
                            //   conformation
};
typedef struct ParameterFit prmset;

#endif
