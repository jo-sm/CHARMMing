#include <time.h>
#include <stdio.h>
#include "Timings.h"

#include "pmeRecipDS.h"

/***=======================================================================***/
/*** InitExecon: initialize a structure to keep timing information during  ***/
/***             an mdgx run.                                              ***/
/***=======================================================================***/
void InitExecon(execon *tm)
{
  tm->bonds = 0.0;
  tm->plist = 0.0;
  tm->nbInt = 0.0;
  tm->nbDirAll = 0.0;
  tm->nbBsp = 0.0;
  tm->nbPtM = 0.0;
  tm->nbFFT = 0.0;
  tm->nbCnv = 0.0;
  tm->nbMtP = 0.0;
  tm->nbMtM = 0.0;
  tm->nbRecAll = 0.0;
  tm->Setup = 0.0;
  tm->Integ = 0.0;
  tm->Write = 0.0;
#ifdef MPI
  tm->mpiMeshPack = 0.0;
  tm->mpiMeshWait = 0.0;
#endif
  gettimeofday(&tm->t0, NULL);
}

/***=======================================================================***/
/*** mdgxStartTimer: record the time just before beginning some piece of   ***/
/***                 code.                                                 ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   tm:     the execon struct in which the timings information is kept  ***/
/***=======================================================================***/
void mdgxStartTimer(execon *tm)
{
  gettimeofday(&tm->tti, NULL);
}

/***=======================================================================***/
/*** mdgxStopTimer: record the time just after finishing some piece of     ***/
/***                code.  Before exiting, this function automatically     ***/
/***                sets the initial time field to the final time field,   ***/
/***                as if mdgxStartTime were implicitly called, so that    ***/
/***                multiple calls to mdgxStopTimer() can be used to       ***/
/***                record a number of successive time intervals without   ***/
/***                redundant calls to mdgxStartTimer().                   ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   tm:     the execon struct in which the timings information is kept  ***/
/***=======================================================================***/
double mdgxStopTimer(execon *tm)
{
  double dt;

  gettimeofday(&tm->ttf, NULL);
  dt = tm->ttf.tv_sec - tm->tti.tv_sec +
    (1.0e-6)*(tm->ttf.tv_usec - tm->tti.tv_usec);
  tm->tti = tm->ttf;

  return dt;
}

/***=======================================================================***/
/*** PrintTimingData: this function prints out all data related to timings ***/
/***                  accumulated over the course of the run.              ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   tm:       the timing data structure                                 ***/
/***   outp:     the output file                                           ***/
/***=======================================================================***/
void PrintTimingData(execon *tm, reccon *rcinp, FILE *outp)
{
  double tt;

  /*** Compute total time spent in this run ***/
  tt = tm->bonds + tm->plist + tm->nbInt + tm->nbBsp + tm->nbPtM + tm->nbFFT +
    tm->nbMtP + tm->nbMtM + tm->Setup + tm->Integ + tm->Write;
#ifdef MPI
  tt += tm->mpiMeshWait + tm->mpiMeshPack;
#endif

  fprintf(outp, " Timings are expressed as a percentage of total tracked "
	  "CPU execution time.\n\n"
	  " Segment                 Time         Percentage \n"
	  " ----------------------  ----------   -----------\n");
  fprintf(outp, " Setup:                  %10.2lf   ( %6.2lf %%)\n", tm->Setup,
	  100.0*tm->Setup/tt);
  fprintf(outp, " Bonded Interactions:    %10.2lf   ( %6.2lf %%)\n", tm->bonds,
	  100.0*tm->bonds/tt);
  fprintf(outp, " Pair List Build:        %10.2lf   ( %6.2lf %%)\n", tm->plist,
	  100.0*tm->plist/tt);
  fprintf(outp, " Nonbonded Interactions: %10.2lf   ( %6.2lf %%)\n", tm->nbInt,
	  100.0*tm->nbInt/tt);
  fprintf(outp, " B-Spline Computation:   %10.2lf   ( %6.2lf %%)\n", tm->nbBsp,
	  100.0*tm->nbBsp/tt);
  fprintf(outp, " Particle -> Mesh:       %10.2lf   ( %6.2lf %%)\n", tm->nbPtM,
	  100.0*tm->nbPtM/tt);
#ifdef MPI
  fprintf(outp, " Mesh Transmission:      %10.2lf   ( %6.2lf %%)\n",
	  tm->mpiMeshWait, 100.0*tm->mpiMeshWait/tt);
  fprintf(outp, " Mesh Packaging:         %10.2lf   ( %6.2lf %%)\n",
	  tm->mpiMeshPack, 100.0*tm->mpiMeshPack/tt);
#endif
  fprintf(outp, " Convolution kernel:     %10.2lf   ( %6.2lf %%)\n", tm->nbCnv,
	  100.0*tm->nbCnv/tt);
  fprintf(outp, " Fast Fourier Transform: %10.2lf   ( %6.2lf %%)\n", tm->nbFFT,
	  100.0*tm->nbFFT/tt);
  if (rcinp->nlev > 1) {
    fprintf(outp, " Mesh -> Mesh:           %10.2lf   ( %6.2lf %%)\n",
	    tm->nbMtM, 100.0*tm->nbMtM/tt);
  }
  fprintf(outp, " Mesh -> Particle:       %10.2lf   ( %6.2lf %%)\n", tm->nbMtP,
	  100.0*tm->nbMtP/tt);
  fprintf(outp, " Integration:            %10.2lf   ( %6.2lf %%)\n", tm->Integ,
	  100.0*tm->Integ/tt);
  fprintf(outp, " Output Printing:        %10.2lf   ( %6.2lf %%)\n", tm->Write,
	  100.0*tm->Write/tt);
  fprintf(outp, " Total CPU Time:         %10.2lf   ( %6.2lf %%)\n", tt,
	  100.0);

  /*** Reset the timers, in case this file is but one ***/
  /*** in a series that is to be printed.             ***/
  InitExecon(tm);
  mdgxStartTimer(tm);
}
