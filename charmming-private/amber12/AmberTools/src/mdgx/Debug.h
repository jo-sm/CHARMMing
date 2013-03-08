#ifndef DebugFunctions
#define DebugFunctions

#include "CrdManipDS.h"
#include "CellManipDS.h"

void PrintCellContents(cellgrid *CG, char* outname, char* varname);

void FindAllInstances(cellgrid *CG, coord *crd, prmtop *tp, int aid);

void CheckCellContents(cellgrid *CG, coord *crd, prmtop *tp, int seekEP,
		       int chkbounds, int chkforces, int StopOnError);

#ifdef MPI
void PrintRecvInfo(int maxsize, char* tname, int recver, int sender, int tag,
		   MPI_Comm comm);

void PrintSendInfo(int maxsize, char* tname, int recver, int sender, int tag,
		   MPI_Comm comm, void* data, int verbose);
#endif

void Torque(prmtop *tp, coord *crd, int resid);

void PrintResidueExclusions(prmtop *tp, int rid);

void PrintResidueForces(cellgrid *CG, coord *crd, prmtop *tp, int rid);

void CellChecksum(cellgrid *CG, coord *crd, int nsctr, int doloc, int dovel,
                  int doploc, int dopvel, char* tagmsg);

#endif
