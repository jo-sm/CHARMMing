#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mdgxVector.h"
#include "Topology.h"
#include "CellManip.h"

#include "CrdManipDS.h"

/***=======================================================================***/
/*** PrintCellContents: loop over all cells and print the contents of the  ***/
/***                    primary sectors to a file.                         ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   CG:       the cell grid                                             ***/
/***   outname:  the name of the output file (MatLab format)               ***/
/***   varname:  the name of the variable for accessing data in MatLab     ***/
/***=======================================================================***/
void PrintCellContents(cellgrid *CG, char* outname, char* varname)
{
  int i, j;
  cell *C;
  FILE *outp;

  outp = fopen(outname, "w");
  const int ncell = CG->ng[0]*CG->ng[1]*CG->ng[2];
  for (i = 0; i < ncell; i++) {
    C = &CG->data[i];
    fprintf(outp, "%% Cell %3d: %4d atoms\n", i, C->nr[0]);
    fprintf(outp, "%s(1:%d,:,%d) = [\n", varname, C->nr[0], i);
    for (j = 0; j < C->nr[0]; j++) {
      fprintf(outp, "%4d %4d %16.10lf %16.10lf %16.10lf %16.10lf\n",
	      C->data[j].id, C->data[j].lj, C->data[j].q, C->data[j].loc[0],
	      C->data[j].loc[1], C->data[j].loc[2]);
    }
    fprintf(outp, "];\n");
  }
  fclose(outp);
}

/***=======================================================================***/
/*** FindAllInstances: find all instances of an atom in all cells.         ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   CG:      the cell grid                                              ***/
/***   crd:     the coordinates                                            ***/
/***   aid:     the atom number to search for                              ***/
/***=======================================================================***/
void FindAllInstances(cellgrid *CG, coord *crd, int aid)
{
  int ci, cj, ck, j, k;
  cell *C;
  atomc *atmt;

  /*** Print some basic information about the atom ***/
  printf("Finding all instances of atom %d:\n", aid);
  printf("Most recent global location: [ %16.10lf %16.10lf %16.10lf ]\n",
	 crd->loc[3*aid], crd->loc[3*aid+1], crd->loc[3*aid+2]);

  /*** Loop over all cells and find instances in any sector ***/
  for (ci = 0; ci < CG->ng[0]; ci++) {
    for (cj = 0; cj < CG->ng[1]; cj++) {
      for (ck = 0; ck < CG->ng[2]; ck++) {
	C = &CG->map[ci][cj][ck];
	for (j = 0; j < 8; j++) {
	  for (k = 0; k < C->nr[j]; k++) {
	    atmt = &C->map[j][k];
	    if (atmt->id == aid) {
	      printf("Atom %6d found in cell [%2d][%2d][%2d]:  ", aid, ci, cj,
		     ck);
	      printf("Sector %d -> Element %4d of %4d (max %4d)\n", j, k,
		     C->nr[j], C->maxatom);
	      printf("Cell Loc: [ %16.10lf %16.10lf %16.10lf ]\n",
		     atmt->loc[0], atmt->loc[1], atmt->loc[2]);
	      printf("Cell Lim: [ %16.10lf %16.10lf %16.10lf ] ->\n"
		     "          [ %16.10lf %16.10lf %16.10lf ]\n",
		     C->orig[0], C->orig[1], C->orig[2],
		     C->orig[0] + CG->celldim[0] + CG->celldim[3],
		     C->orig[1] + CG->celldim[1] + CG->celldim[3],
		     C->orig[2] + CG->celldim[2] + CG->celldim[3]);
	    }
	  }
	}
      }
    }
  }
}

/***=======================================================================***/
/*** PrimeCCCErrMsg: initialize the error message for CheckCellContents()  ***/
/***                 below.  Critical information supplied includes the    ***/
/***                 process, system, and atom identifications.            ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   CG:           the cell grid                                         ***/
/***   C:            the cell within the cell grid which has the error     ***/
/***   sctr:         the sector of the cell wherein the problem is found   ***/
/***   cidx:         the index of the troublesome atom within the cell     ***/
/***                 (the local variable tpidx gives the "absolute" index  ***/
/***                 within the actual topology)                           ***/
/***   tp:           the topology                                          ***/
/***   errname:      the name of the error                                 ***/
/***=======================================================================***/
static int PrimeCCCErrMsg(cellgrid *CG, cell *C, int sctr, int cidx,
			  prmtop *tp, char* errname, char* ErrMsg)
{
  int tpidx, residx;
  char *errtmp;

  tpidx = (sctr >= 0) ? C->map[sctr][cidx].id : cidx;
  residx = LocateResID(tp, tpidx, 0, tp->nres);
  sprintf(ErrMsg, ">> %s on process %4d of system %3d:\n", errname, CG->tid,
	  CG->sysID);
  errtmp = &ErrMsg[strlen(ErrMsg)];
  sprintf(errtmp, ">> Atom %6d (%.4s %4d %.4s)", tpidx,
	  &tp->ResNames[4*residx], residx, &tp->AtomNames[4*tpidx]);
  if (sctr >= 0) {
    errtmp = &ErrMsg[strlen(ErrMsg)];
    sprintf(errtmp, " in cell [ %3d %3d %3d ][ %d ] -> %4d:", C->gbin[0],
	    C->gbin[1], C->gbin[2], sctr, cidx);
  }
  errtmp = &ErrMsg[strlen(ErrMsg)];
  sprintf(errtmp, "\n");

  return strlen(ErrMsg);
}

/***=======================================================================***/
/*** CheckCellContents: this debugging function performs a complex assert, ***/
/***                    verifying that the contents of each cell sector    ***/
/***                    really belong there.  This function will work in   ***/
/***                    parallel implementations, and should serve as a    ***/
/***                    useful debugging tool.                             ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   CG:           the cell grid                                         ***/
/***   crd:          the coordinates                                       ***/
/***   tp:           the topology                                          ***/
/***   seekEP:       flag to activate seeking of extra points in certain   ***/
/***                 tests                                                 ***/
/***   chkbounds:    flag to call for bounds checking on all cells--set to ***/
/***                 zero to turn this off if constraints have just been   ***/
/***                 applied or Verlet integration has just occured and    ***/
/***                 the UpdateCells() function of the CellManip library   ***/
/***                 has not yet been called                               ***/
/***   chkforces:    flag to call for force checking, which will report    ***/
/***                 any unusually large forces                            ***/
/***   StopOnError:  flag to have the program abort if errors are detected ***/
/***=======================================================================***/
void CheckCellContents(cellgrid *CG, coord *crd, prmtop *tp, int seekEP,
		       int chkbounds, int chkforces, int StopOnError)
{
  int i, j, k, m, cp3, problem, halt, msgidx, nduploc;
  int* cellplace;
  int* buffplace;
  int* duploc;
  double fmag;
  double climit[3];
  char* ErrMsg;
  char *errtmp;
  atomc *catm;
  cell *C, *CNew;

  /*** Allocate memory for the error message, ***/
  /*** for duplicate atoms thus far located,  ***/
  /*** and other needs                        ***/
  ErrMsg = (char*)malloc(1024*sizeof(char));
  duploc = (int*)malloc(CG->maxatom*sizeof(int));

  /*** Check to see that no atoms are outside of their cell boundaries ***/
  halt = 0;
  if (chkbounds == 1) {
    for (i = 0; i < CG->MyCellCount; i++) {
      C = &CG->data[CG->MyCellDomain[i]];

      /*** Shortcut for cell limits ***/
      for (j = 0; j < 3; j++) {
	climit[j] = C->orig[j] + CG->celldim[j] + CG->celldim[3];
      }
      for (j = 0; j < 8; j++) {
	catm = C->map[j];
	for (k = 0; k < C->nr[j]; k++) {

	  /*** Determine if the atom is outside the cell boundary ***/
	  problem = 0;
	  for (m = 0; m < 3; m++) {
	    if (catm[k].loc[m] < C->orig[m] || catm[k].loc[m] >= climit[m]) {
	      problem = 1;
	    }
	  }
	  if (problem == 1) {
	    msgidx = PrimeCCCErrMsg(CG, C, j, k, tp, "Bounds alert", ErrMsg);
	    errtmp = &ErrMsg[msgidx];
	    sprintf(errtmp, "   Atom location: %10.5lf %10.5lf %10.5lf\n",
		    catm[k].loc[0], catm[k].loc[1], catm[k].loc[2]);
	    errtmp = &ErrMsg[strlen(ErrMsg)];
	    sprintf(errtmp, "   Cell limits: [ %10.5lf %10.5lf %10.5lf ] x\n"
		            "                [ %10.5lf %10.5lf %10.5lf ]\n\n",
		    C->orig[0], C->orig[1], C->orig[2], climit[0], climit[1],
		    climit[2]);
	    printf("%s\n", ErrMsg);
	    halt = 1;
	  }
	}
      }
    }
  }

  /*** Check for duplicate atoms in any cells ***/
  nduploc = 0;
  for (i = 0; i < CG->MyCellCount; i++) {
    C = &CG->data[CG->MyCellDomain[i]];
    for (j = 0; j < 8; j++) {
      catm = C->map[j];
      for (k = 0; k < C->nr[j]-1; k++) {
	if (DetectInIVec(catm[k].id, duploc, nduploc) >= 0) {
	  continue;
	}
	problem = 0;
	for (m = k+1; m < C->nr[j]; m++) {
	  if (catm[k].id == catm[m].id) {
	    if (problem == 0) {
	      msgidx = PrimeCCCErrMsg(CG, C, 0, k, tp, "Duplicate atom",
				      ErrMsg);
	      duploc[nduploc] = catm[k].id;
	      nduploc++;
	    }
	    errtmp = &ErrMsg[msgidx];
	    sprintf(errtmp, "   Also occupies index %3d.\n", m);
	    problem = 1;
	    halt = 1;
	  }
	}
	if (problem == 1) {
	  printf("%s\n", ErrMsg);
	}
      }
    }
  }

  /*** Check to see that there are no unusually large forces ***/
  if (chkforces == 1) {
    for (i = 0; i < CG->MyCellCount; i++) {
      C = &CG->data[CG->MyCellDomain[i]];
      for (j = 0; j < C->nr[0]; j++) {
	fmag = C->data[j].frc[0]*C->data[j].frc[0];
	fmag += C->data[j].frc[1]*C->data[j].frc[1];
	fmag += C->data[j].frc[2]*C->data[j].frc[2];
	fmag = sqrt(fmag);
	if (fmag > 10000.0) {
	  k = LocateResID(tp, C->data[j].id, 0, tp->nres);
	  msgidx = PrimeCCCErrMsg(CG, C, 0, k, tp, "Large force", ErrMsg);
	  errtmp = &ErrMsg[msgidx];
	  sprintf(errtmp, "   Force = [ %16.4lf %16.4lf %16.4lf ]\n",
		  C->data[j].loc[0], C->data[j].loc[1], C->data[j].loc[2]);
	  printf("%s\n", ErrMsg);
	  halt = 1;
	}
      }
    }
  }

  /*** Check to see that all atoms are accounted for ***/
  cellplace = (int*)malloc(crd->natom*sizeof(int));
  SetIVec(cellplace, crd->natom, -1);
  for (i = 0; i < CG->MyCellCount; i++) {
    C = &CG->data[CG->MyCellDomain[i]];
    catm = C->data;
    for (m = 0; m < C->nr[0]; m++) {

      /*** First check if this atom has already been located ***/
      if (catm[m].id < 0 || catm[m].id >= crd->natom) {
	printf("Error: atom [ %2d %2d %2d ]->[ 0, %4d ] ID = %d!\n",
	       i, j, k, m, catm[m].id);
	exit(1);
      }
      cp3 = catm[m].id;
      if (cellplace[cp3] >= 0) {
	CNew = &CG->data[cellplace[cp3]];
        msgidx = PrimeCCCErrMsg(CG, C, 0, k, tp, "Multiple cell locations",
				ErrMsg);
        errtmp = &ErrMsg[msgidx];
	sprintf(errtmp, "   Also located in [ %3d %3d %3d ]!\n", CNew->gbin[0],
		CNew->gbin[1], CNew->gbin[2]);
	printf("%s\n", ErrMsg);
	halt = 1;
      }

      /*** Mark this atom as located ***/
      cellplace[cp3] = C->gbin[3];
    }
  }

  /*** Check to see that all cells contain properly ordered lists ***/
  for (i = 0; i < CG->MyCellCount; i++) {
    C = &CG->data[CG->MyCellDomain[i]];
    for (j = 0; j < 8; j++) {
      catm = C->map[j];
      for (k = 1; k < C->nr[j]; k++) {
	if (catm[k].id < catm[k-1].id) {
	  PrimeCCCErrMsg(CG, C, j, k, tp, "Improper order", ErrMsg);
	  printf("%s", ErrMsg);
	  PrimeCCCErrMsg(CG, C, j, k-1, tp, "   Out of order relative to",
			 ErrMsg);
	  printf("%s\n", ErrMsg);
	  halt = 1;
	}
      }
    }
  }

  /*** Check to see that all atoms are accounted for ***/
  buffplace = (int*)malloc(crd->natom*sizeof(int));
#ifdef MPI
  MPI_Reduce(cellplace, buffplace, crd->natom, MPI_INT, MPI_MAX, 0,
	     CG->dspcomm);
#endif
  ReflectIVec(cellplace, buffplace, crd->natom);
  if (CG->tid == 0) {
    for (i = 0; i < tp->natom; i++) {
      if (cellplace[i] == -1 && (tp->Masses[i] > 1.0e-8 || seekEP == 1)) {
	PrimeCCCErrMsg(CG, C, -1, i, tp, "Nonassigned atom", ErrMsg);
	printf("%s", ErrMsg);
      }
    }
  }
  if (halt == 1 && StopOnError == 1) {
    printf("Exiting due to previous errors.\n");
    exit(1);
  }

#ifdef MPI
  /*** Sync processes ***/
  MPI_Barrier(CG->dspcomm);
#endif

  /*** Free allocated memory ***/
  free(cellplace);
  free(buffplace);
  free(ErrMsg);
  free(duploc);
}

#ifdef MPI
/***=======================================================================***/
/*** PrintRecvInfo: print the information relating to a recv posting.      ***/
/***                This function is intended to speed the debugging       ***/
/***                process by encapsulating this sort of report into a    ***/
/***                form that can be added conveniently anywhere in the    ***/
/***                main code.                                             ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   maxsize:     the maximum size of the incoming message               ***/
/***   tname:       string defining the data type (for example, "INT")     ***/
/***   recver:      the receiving process                                  ***/
/***   sender:      the process from which to expect the message           ***/
/***   tag:         the message tag                                        ***/
/***   comm:        the communicator over which the recv passes            ***/
/***=======================================================================***/
void PrintRecvInfo(int maxsize, char* tname, int recver, int sender, int tag,
		   MPI_Comm comm)
{
  printf("Posting recv %4d -> %4d (%10d): Max  %6d of %s over %10d\n", sender,
	 recver, tag, maxsize, tname, (int)comm);
}

/***=======================================================================***/
/*** PrintSendInfo: print the information relating to a send posting.      ***/
/***                This function is intended to speed the debugging       ***/
/***                process by encapsulating this sort of report into a    ***/
/***                form that can be added conveniently anywhere in the    ***/
/***                main code.                                             ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   actsize:     the actual size of the incoming message                ***/
/***   tname:       string defining the data type (for example, "INT")     ***/
/***   recver:      the receiving process                                  ***/
/***   sender:      the process from which to expect the message           ***/
/***   tag:         the message tag                                        ***/
/***   comm:        the communicator over which the recv passes            ***/
/***   data:        buffer containing the message                          ***/
/***   verbose:     flag to activate printing of ALL message data          ***/
/***=======================================================================***/
void PrintSendInfo(int actsize, char* tname, int recver, int sender, int tag,
		   MPI_Comm comm, void* data, int verbose)
{
  int i, outsize, typID;
  char *ctmp;
  char* output;

  /*** Considerations for printing ATOMB type data ***/
  outsize = 128;
  if (strcmp(tname, "ATOMB") == 0) {
    outsize += actsize*48;
    typID = 1;
  }
  else if (strcmp(tname, "ATOMBV") == 0) {
    outsize += actsize*80;
    typID = 2;
  }
  else if (strcmp(tname, "ATOMBX") == 0) {
    outsize += actsize*220;
    typID = 3;
  }
  else {

    /*** Verbose printing is not supported for this data type ***/
    verbose = 0;
  }
  output = (char*)malloc(outsize*sizeof(char));
  sprintf(output, "Posting send %4d -> %4d (%10d): Size %6d of %s over %10d\n",
	  sender, recver, tag, actsize, tname, (int)comm);
  if (verbose == 0) {

    /*** Print message and free the memory before returning ***/
    printf("%s", output);
    free(output);

    return;
  }

  /*** If we're still here, print the entire send ***/
  ctmp = &output[strlen(output)];
  if (typID == 1) {
    atomb *rdata;
    rdata = (atomb*)data;

    /*** It is assumed that, if the message is one of these ***/
    /*** custom-defined atom-related types, that the first  ***/
    /*** element is not occupied with actual data.          ***/
    for (i = 1; i < actsize+1; i++) {
      sprintf(ctmp, "  %6d %4d %10.4lf %10.4lf %10.4lf\n", rdata[i].id,
	      rdata[i].dreg, rdata[i].loc[0], rdata[i].loc[1],
	      rdata[i].loc[2]);
      ctmp = &output[strlen(output)];
    }
  }
  if (typID == 2) {
    atombv *rdata;
    rdata = (atombv*)data;

    /*** It is assumed that, if the message is one of these ***/
    /*** custom-defined atom-related types, that the first  ***/
    /*** element is not occupied with actual data.          ***/
    for (i = 1; i < actsize+1; i++) {
      sprintf(ctmp, "  %6d %4d %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf"
	      "\n", rdata[i].id, rdata[i].dreg, rdata[i].loc[0],
	      rdata[i].loc[1], rdata[i].loc[2], rdata[i].vel[0],
	      rdata[i].vel[1], rdata[i].vel[2]);
      ctmp = &output[strlen(output)];
    }
  }
  printf("%s", output);

  /*** Free allocated memory ***/
  free(output);
}
#endif

/***=======================================================================***/
/*** Torque: compute the torque for a residue.                             ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   tp:     the topology structure                                      ***/
/***   crd:    the system coordinates                                      ***/
/***   resid:  the residue number                                          ***/
/***=======================================================================***/
static void Torque(prmtop *tp, coord *crd, int resid)
{
  int i, j;
  double com[3], r[3], rxf[3], trq[3];

  /*** Compute center of mass ***/
  com[0] = com[1] = com[2] = 0.0;
  for (i = tp->ResLims[resid]; i < tp->ResLims[resid+1]; i++) {
    for (j = 0; j < 3; j++) {
      com[j] += tp->Masses[i]*crd->loc[3*i+j];
    }
  }

  /*** Compute torque ***/
  trq[2] = trq[1] = trq[0] = 0.0;
  for (i = tp->ResLims[resid]; i < tp->ResLims[resid+1]; i++) {
    for (j = 0; j < 3; j++) {
      r[j] = crd->loc[3*i+j] - com[j];
    }
    CrossP(r, &crd->frc[3*i], rxf);
    for (j = 0; j < 3; j++) {
      trq[j] += rxf[j];
    }
  }

  printf("Torque = [ %12.8lf %12.8lf %12.8lf ]\n", trq[0], trq[1], trq[2]);
}

/***=======================================================================***/
/*** WriteElimPair: write an elimination pair.  This condenses the code in ***/
/***                PrintResidueEliminations below.                        ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   tp:        the topology                                             ***/
/***   A:         integer 1, 2, 3, or 4                                    ***/
/***   nA:        the number of 1:A eliminations                           ***/
/***   listA:     the list of 1:A eliminations                             ***/
/***=======================================================================***/
static void WriteElimPair(prmtop *tp, int A, int nA, nixpr* listA)
{
  int i, xres, yres;

  for (i = 0; i < nA; i++) {
    xres = LocateResID(tp, listA[i].atmX, 0, tp->nres);
    yres = LocateResID(tp, listA[i].atmY, 0, tp->nres);
    printf("  1:%d -> %6d / %6d ( %.4s %.4s / %.4s %.4s )\n", A, listA[i].atmX,
	   listA[i].atmY, &tp->AtomNames[4*listA[i].atmX],
	   &tp->ResNames[4*xres], &tp->AtomNames[4*listA[i].atmY],
	   &tp->ResNames[4*yres]);
  }
}

/***=======================================================================***/
/*** PrintResidueExclusions: print the exclusion list, atom by atom, for   ***/
/***                         the specified residue.                        ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   tp:        the topology                                             ***/
/***   rid:       the residue ID                                           ***/
/***=======================================================================***/
void PrintResidueExclusions(prmtop *tp, int rid)
{
  int i, j, llim, hlim, jatm, jres;

  llim = tp->ResLims[rid];
  hlim = tp->ResLims[rid+1];
  for (i = llim; i < hlim; i++) {

    /*** Print exclusions ***/
    printf("Exclusions for atom %6d ( %.4s %.4s ):\n", i, &tp->AtomNames[4*i],
	   &tp->ResNames[4*rid]);
    for (j = tp->ConExcl[i]; j < tp->ConExcl[i+1]; j++) {
      jatm = tp->ExclList[j];

      /*** Bail out if there are no exclusions ***/
      if (jatm == -1) {
	continue;
      }
      jres = LocateResID(tp, jatm, 0, tp->nres);
      printf("  %6d %.4s %.4s\n", jatm, &tp->AtomNames[4*jatm],
	     &tp->ResNames[4*jres]);
    }

    /*** Print eliminations ***/
    if (tp->EPInserted == 1) {
      printf("Eliminations for atom %6d ( %.4s %.4s ):\n", i,
	     &tp->AtomNames[4*i], &tp->ResNames[4*rid]);
      WriteElimPair(tp, 1, tp->ElimPair[i].n11, tp->ElimPair[i].list11);
      WriteElimPair(tp, 2, tp->ElimPair[i].n12, tp->ElimPair[i].list12);
      WriteElimPair(tp, 3, tp->ElimPair[i].n13, tp->ElimPair[i].list13);
      WriteElimPair(tp, 4, tp->ElimPair[i].n14, tp->ElimPair[i].list14);
      printf("\n");
    }
  }
}

/***=======================================================================***/
/*** PrintResidueForces: print forces on all atoms of a particular residue ***/
/***                     to stdout.  This function searches the cell grid  ***/
/***                     primary sectors for all atoms.                    ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   CG:        the cell grid                                            ***/
/***   tp:        the topology                                             ***/
/***   rid:       the residue ID                                           ***/
/***=======================================================================***/
void PrintResidueForces(cellgrid *CG, coord *crd, prmtop *tp, int rid)
{
  int h, i, j, k, atmcid;
  double *ftmp;

  printf("%% Residue %d:\n", rid);
  printf("%%  Force X       Force Y       Force Z   "
	 "  Velocity X    Velocity Y    Velocity Z   Cell Index "
	 "  Position X    Position Y    Position Z    Charge\n"
	 "------------- ------------- ------------- "
	 "------------- ------------- ------------- --- --- --- "
	 "------------- ------------- ------------- ----------\n");
  for (h = tp->ResLims[rid]; h < tp->ResLims[rid+1]; h++) {
    for (i = 0; i < CG->ng[0]; i++) {
      for (j = 0; j < CG->ng[1]; j++) {
	for (k = 0; k < CG->ng[2]; k++) {
	  atmcid = FindAtomInSector(&CG->map[i][j][k], h, 0, 0);
	  if (atmcid >= 0) {
	    ftmp = CG->map[i][j][k].data[atmcid].frc;
	    printf("%13.7lf %13.7lf %13.7lf %13.7lf %13.7lf %13.7lf", ftmp[0],
		   ftmp[1], ftmp[2], crd->vel[3*h], crd->vel[3*h+1],
		   crd->vel[3*h+2]);
	    ftmp = CG->map[i][j][k].data[atmcid].loc;
	    printf(" %3d %3d %3d %13.7lf %13.7lf %13.7lf %10.7lf\n", i, j, k,
		   ftmp[0], ftmp[1], ftmp[2], CG->map[i][j][k].data[atmcid].q);
	    i = CG->ng[0];
	    j = CG->ng[1];
	    k = CG->ng[2];
	  }
	}
      }
    }
  }
  printf("\n");
}

/***=======================================================================***/
/*** CellChecksum: perform a checksum operation on all cells, over various ***/
/***               values depending on the context.                        ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   CG:        the cell grid                                            ***/
/***   crd:       the coordinates struct                                   ***/
/***   nsctr:     the number of sectors for which to print results (if     ***/
/***              this is set to 1, only the primary sector results will   ***/
/***              be displayed)                                            ***/
/***   do[p]loc:  flag to activate (previous) atom location sums           ***/
/***   do[p]vel:  flag to activate (previous) atom velocity sums           ***/
/***   tagmsg:    message tag for this checksum operation                  ***/
/***=======================================================================***/
void CellChecksum(cellgrid *CG, coord *crd, int nsctr, int doloc, int dovel,
		  int doploc, int dopvel, char* tagmsg)
{
  int i, j, k, atmid, olen;
  double chksum[4];
  char *ctmp;
  char* outline;
  cell *C;

  /*** Allocate memory for output message ***/
  outline = (char*)malloc(32 + (doloc + dovel + doploc + dopvel) *
			  nsctr*21*sizeof(char));

#ifdef MPI
  /*** Sync up all processes ***/
  MPI_Barrier(CG->dspcomm);
#endif
  if (CG->tid == 0) {
    printf("%s\n", tagmsg);
  }
  for (i = 0; i < CG->ncell; i++) {
#ifdef MPI
    MPI_Barrier(CG->dspcomm);
#endif
    if (CG->data[i].CGRank != CG->tid) {
      continue;
    }
    C = &CG->data[i];

    /*** Print the checksums for all sectors ***/
    sprintf(outline, "%4d = [ ", C->gbin[3]);
    olen = strlen(outline);
    ctmp = &outline[olen];
    for (j = 0; j < nsctr; j++) {

      /*** Accumulate checksums ***/
      for (k = 0; k < 4; k++) {
	chksum[k] = 0.0;
      }
      for (k = 0; k < C->nr[j]; k++) {
	atmid = C->map[j][k].id;

	/*** Current coordinates ***/
	if (doloc) {
	  chksum[0] += C->map[j][k].loc[doloc-1];
	}

	/*** Current velocities ***/
	if (dovel) {
	  chksum[1] += crd->vel[3*atmid+dovel-1];
	}

	/*** Previous coordinates ***/
	if (doploc) {
	  chksum[2] += crd->prvloc[3*atmid+doploc-1];
	}

	/*** Previous velocities ***/
	if (dopvel) {
	  chksum[3] += crd->prvvel[3*atmid+doploc-1];
	}
      }

      /*** Print the checksums ***/
      if (doloc) {
	sprintf(ctmp, "%20.12e ", chksum[0]);
	olen = strlen(outline);
	ctmp = &outline[olen];
      }
      if (dovel) {
	sprintf(ctmp, "%20.12e ", chksum[1]);
	olen = strlen(outline);
	ctmp = &outline[olen];
      }
      if (doploc) {
	sprintf(ctmp, "%20.12e ", chksum[2]);
	olen = strlen(outline);
	ctmp = &outline[olen];
      }
      if (dopvel) {
	sprintf(ctmp, "%20.12e ", chksum[3]);
	olen = strlen(outline);
	ctmp = &outline[olen];
      }
    }
    sprintf(ctmp, "];\n");
    printf("%s", outline);
  }
}
