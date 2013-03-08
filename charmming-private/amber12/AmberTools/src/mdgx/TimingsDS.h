#ifndef TimingsStructs
#define TimingsStructs

#include <sys/time.h>

struct ExecutionControl {
  struct timeval t0;
  struct timeval tC;
  struct timeval tti;
  struct timeval ttf;
  double bonds;
  double plist;
  double nbInt;
  double nbDirAll;
  double nbBsp;
  double nbPtM;
  double nbFFT;
  double nbCnv;
  double nbMtP;
  double nbMtM;
  double nbRecAll;
  double Setup;
  double Integ;
  double Write;
#ifdef MPI
  double mpiMeshWait;
  double mpiMeshPack;
#endif
};
typedef struct ExecutionControl execon;

#endif
