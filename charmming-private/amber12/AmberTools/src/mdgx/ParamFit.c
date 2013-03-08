#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mdgxVector.h"
#include "Matrix.h"
#include "ParamFit.h"
#include "Topology.h"
#include "CrdManip.h"
#include "Parse.h"
#include "ChargeFit.h"
#include "VirtualSites.h"
#include "Manual.h"

#include "TrajectoryDS.h"

/***=======================================================================***/
/*** AssignChargeGroups: this function computes the charge groups to be    ***/
/***                     manipulated during parameter optimization.        ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:   parameter optimization controls                               ***/
/***=======================================================================***/
static void AssignChargeGroups(prmset *mp, prmtop *tp)
{
  int i, j, natm;
  int* atmmask;
  double qsum;
  coord crd;

  /*** Find the atoms corresponding to each variable charge group ***/
  crd = CreateCoord(tp->natom);
  for (i = 0; i < mp->nqvar; i++) {
    atmmask = ParseAmbMask(mp->qvar[i].maskstr, tp, &crd);
    natm = ISum(atmmask, tp->natom);
    mp->qvar[i].natm = natm;
    if (natm > 0) {
      mp->qvar[i].atoms = (int*)malloc(natm*sizeof(int));
      natm = 0;
      for (j = 0; j < tp->natom; j++) {
	if (atmmask[j] == 1) {
	  mp->qvar[i].atoms[natm] = j;
	  natm++;
	}
      }
    }
    free(atmmask);
  }
  DestroyCoord(&crd);

  /*** If there is only one group of variable charges, ***/ 
  /*** charge optimization is merely about fulfilling  ***/
  /*** the total charge constraint.                    ***/
  if (mp->nqvar == 1) {
    for (i = 0; i < mp->qvar[0].natm; i++) {
      tp->Charges[mp->qvar[0].atoms[i]] = 0.0;
    }
    qsum = DSum(tp->Charges, tp->natom)/mp->qvar[0].natm;
    for (i = 0; i < mp->qvar[0].natm; i++) {
      tp->Charges[mp->qvar[0].atoms[i]] = -qsum;
    }
    mp->nqvar = 0;

    return;
  }
}

/***=======================================================================***/
/*** IntrExcluded: returns 0 if there is an exclusion (1:1, 1:2, 1:3, or   ***/
/***               1:4) that prevents a full nonbonded interaction from    ***/
/***               being computed, 1 otherwise.                            ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   aatm:   the first atom of the pair                                  ***/
/***   batm:   the second atom of the pair                                 ***/
/***   tp:     the energy function                                         ***/
/***=======================================================================***/
static int IntrExcluded(int aatm, int batm, prmtop *tp)
{
  int i;

  /*** Double-check the exclusion list ***/
  for (i = tp->ConExcl[aatm]; i < tp->ConExcl[aatm+1]; i++) {
    if (tp->ExclList[i] == batm) {
      return 0;
    }
  }
  for (i = tp->ConExcl[batm]; i < tp->ConExcl[batm+1]; i++) {
    if (tp->ExclList[i] == aatm) {
      return 0;
    }
  }

  /*** Check to see that there are vdW or electrostatic interactions ***/
  if ((tp->LJIdx[aatm] >= 0 && tp->LJIdx[batm] >= 0) ||
      (fabs(tp->Charges[aatm]) >= 1.0e-8 &&
       fabs(tp->Charges[batm]) >= 1.0e-8)) {
    return 1;
  }
  else {
    return 0;
  }
}

/***=======================================================================***/
/*** CreateNBList: this function makes a pairlist under the assumptions    ***/
/***               that there are no cutoffs.  The pairlist applies to     ***/
/***               calculations with this library's simplified bond,       ***/
/***               angle, and dihedral calculations.                       ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   tp:     the energy function                                         ***/
/***=======================================================================***/
static imat CreateNBList(prmtop *tp)
{
  int i, j;
  int *itmp;
  imat plist;

  /*** Allocate the pair list ***/
  plist = CreateImat(tp->natom, tp->natom+1);
  for (i = 0; i < tp->natom-1; i++) {

    /*** The first index is the number of pairs for this atom ***/
    itmp = plist.map[i];
    itmp[0] = 0;
    for (j = i+1; j < tp->natom; j++) {
      itmp[itmp[0]+1] = j;
      itmp[0] += IntrExcluded(i, j, tp);
    }
  }

  return plist;
}

/***=======================================================================***/
/*** ConfigNrg: test the energies of a variety of configurations given an  ***/
/***            energy function (in the form of a prmtop struct) and a     ***/
/***            trajectory of the conformations (in the form of a dmat     ***/
/***            matrix with one conformation per row).  All energies are   ***/
/***            computed without cutoffs for conformations that are        ***/
/***            assumed to exist isolated in space.                        ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   tp:     the energy function                                         ***/
/***   plist:  the pair list for all configurations                        ***/
/***   crds:   coordinates of all conformations                            ***/
/***   nrg:    the energies of all conformations (returned)                ***/
/***=======================================================================***/
static void ConfigNrg(prmtop *tp, imat *plist, dmat *crds, double *nrg)
{
  int i, j, k, npair, katm;
  int *itmp;
  double dx, dy, dz, invr, invr6, atmx, atmy, atmz, atmq;
  double *loctmp, *ljutmp;

  /*** Loop over all conformations and compute energies ***/
  for (i = 0; i < crds->row; i++) {

    /*** Initialize energy ***/
    nrg[i] = 0.0;

    /*** Compute all nonbonded pairs ***/
    loctmp = crds->map[i];
    for (j = 0; j < tp->natom; j++) {
      npair = plist->map[j][0];
      itmp = &plist->map[j][1];
      atmx = loctmp[3*j];
      atmy = loctmp[3*j+1];
      atmz = loctmp[3*j+2];
      atmq = BIOQ*tp->Charges[j];
      if (tp->LJIdx[j] >= 0) {
	ljutmp = tp->LJutab.map[tp->LJIdx[j]];
      }
      for (k = 0; k < npair; k++) {

	/*** Distance ***/
	katm = itmp[k];
	dx = loctmp[3*katm] - atmx;
	dy = loctmp[3*katm+1] - atmy;
	dz = loctmp[3*katm+2] - atmz;
	invr = 1.0/sqrt(dx*dx + dy*dy + dz*dz);

	/*** Electrostatics ***/
	nrg[i] += invr*atmq*tp->Charges[katm];

	/*** van der Waals ***/
	if (tp->LJIdx[j] >= 0 && tp->LJIdx[katm] >= 0) {
	  invr6 = invr*invr*invr;
	  invr6 *= invr6;
	  nrg[i] += (ljutmp[2*tp->LJIdx[katm]]*invr6 -
		     ljutmp[2*tp->LJIdx[katm]+1])*invr6;
	}
      }
    }
  }
}

/***=======================================================================***/
/*** ReadSystemConf: read a set of system conformations and store them in  ***/
/***                 a matrix.                                             ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp: fitting control data (contains the coordinate file name and     ***/
/***       records the number of fitting points)                           ***/
/***=======================================================================***/
static dmat ReadSystemConf(prmset *mp, prmtop *tp)
{
  int i, j;
  double *loctmp, *aloc, *bloc, *cloc, *dloc, *eploc, *epgloc, *agloc;
  FILE *inp;
  dmat crds;
  expt *tmr;

  /*** Open the file of conformations ***/
  if ((inp = fopen(mp->sysconf, "r")) == NULL) {
    printf("ReadSystemConf >> Error.  System file %s not found.\n",
	   mp->sysconf);
    exit(1);
  }

  /*** Read the number of points ***/
  fscanf(inp, "%d%d", &i, &mp->nconf);
  if (i != tp->norigatom) {
    printf("ReadSystemConf >> Error.  Topology contains %d atoms, system\n"
	   "ReadSystemConf >> conformations have %d atoms.\n", tp->norigatom,
	   i);
    exit(1);
  }

  /*** Allocate memory for all conformations ***/
  crds = CreateDmat(mp->nconf, 3*tp->natom, 0);

  /*** Read the conformations ***/
  for (i = 0; i < mp->nconf; i++) {
    loctmp = crds.map[i];
    for (j = 0; j < 3*tp->norigatom; j++) {
      fscanf(inp, "%lf", &loctmp[j]);
    }
    if (tp->EPInserted == 1) {
      for (j = tp->natom-1; j >= 0; j--) {
	if (tp->OldAtomNum[j] >= 0) {
	  loctmp[3*j] = loctmp[3*tp->OldAtomNum[j]];
	  loctmp[3*j+1] = loctmp[3*tp->OldAtomNum[j]+1];
	  loctmp[3*j+2] = loctmp[3*tp->OldAtomNum[j]+2];
	}
	else {
          loctmp[3*j] = 0.0;
          loctmp[3*j+1] = 0.0;
          loctmp[3*j+2] = 0.0;
	}
      }
    }
  }

  /*** Place extra points on each conformation ***/
  for (i = 0; i < tp->nxtrapt; i++) {
    tmr = &tp->xtrapts[i];
    for (j = 0; j < mp->nconf; j++) {
      loctmp = crds.map[j];
      aloc = &loctmp[3*tmr->fr1];
      agloc = aloc;
      bloc = &loctmp[3*tmr->fr2];
      if (tmr->frstyle > 1) {
	cloc = &loctmp[3*tmr->fr3];
      }
      if (tmr->frstyle == 6) {
	dloc = &loctmp[3*tmr->fr4];
      }
      eploc = &loctmp[3*tmr->atomid];
      epgloc = eploc;
      XptLocator(aloc, agloc, bloc, cloc, dloc, eploc, epgloc, tmr);
    }
  }

  return crds;
}

/***=======================================================================***/
/*** ReadTargetVal: read a list of target numbers, one per conformation.   ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:     fitting control data (contains the file name and the number ***/
/***           of fitting points)                                          ***/
/***   fname:  the name of the file to read                                ***/
/***=======================================================================***/
static double* ReadTargetVal(prmset *mp, char* fname)
{
  int i, npt;
  double* etrg;
  char line[MAXLINE];
  FILE *inp;
  cmat lwords;

  /*** Allocate ***/
  etrg = (double*)malloc(mp->nconf*sizeof(double));

  /*** Read energies ***/
  if ((inp = fopen(fname, "r")) == NULL) {
    printf("ReadTargetVal >> Error.  Energy file %s not found.\n", fname);
    exit(1);
  }
  npt = 0;
  etrg = (double*)malloc(mp->nconf*sizeof(double));
  while (fgets(line, MAXLINE, inp) != NULL) {
    RemoveWhiteSpace(line, MAXLINE);
    if ((line[0] >= '0' && line[0] <= '9') ||
	line[0] == '.' || line[0] == '-') {
      RemoveComments(line);
      lwords = ParseWords(line);
      for (i = 0; i < lwords.row; i++) {
	etrg[npt] = atof(lwords.map[i]);
	npt++;
      }
    }
  }
  if (npt != mp->nconf) {
    printf("ReadTargetVal >> Error.  %d conformations specified, but only "
	   "%d\nReadTargetVal >> conformation energies found.\n",
	   mp->nconf, npt);
    exit(1);
  }

  return etrg;
}

/***=======================================================================***/
/*** ComputeNrgWeights: compute weights for each conformation based on a   ***/
/***                    Boltzmann distribution.  More complex weighting    ***/
/***                    schemes will be made avaiable by reading a second  ***/
/***                    file of numbers for each conformation.             ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:   parameter fitting information                                 ***/
/***   etrg:      vector of energies for each conformation                 ***/
/***=======================================================================***/
static double* ComputeNrgWeights(prmset *mp, double* etrg)
{
  int i;
  double tmass;
  double* ewt;

  /*** Check for a user-specified weight scheme ***/
  if (mp->syswt[0] != '\0') {
    ewt = ReadTargetVal(mp, mp->syswt);
    return ewt;
  }

  /*** Compute weights by standard Boltzmann distribution ***/
  ewt = (double*)malloc(mp->nconf*sizeof(double));
  tmass = 0.0;
  for (i = 0; i < mp->nconf; i++) {
    ewt[i] = exp(-etrg[i]/(GASCNST*mp->wttemp));
    tmass += ewt[i];
  }
  tmass = mp->nconf/tmass;
  for (i = 0; i < mp->nconf; i++) {
    ewt[i] *= tmass;
  }

  return ewt;
}

/***=======================================================================***/
/*** SeedNewCharges: changes the charges in a topology in accord with the  ***/
/***                 variable charge groups and seeds found in the         ***/
/***                 parameter input.                                      ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/

/***   mp:   parameter fitting information                                 ***/
/***=======================================================================***/
static int SeedNewCharges(prmtop *tp, prmtop *newtp, prmset *mp)
{
  int i, j, pivot;
  int *atmlist;
  double qsum;
  double* qval;

  /*** Figure out the charges on each group ***/
  qval = (double*)malloc(mp->nqvar*sizeof(double));
  for (i = 0; i < mp->nqvar; i++) {
    qval[i] = mp->qvar[i].target -
      (0.5*(mp->nqseed-1) - mp->qseed[i])*mp->qspc;
  }

  /*** The last group of variable charges is used ***/
  /*** to balance the sum of previous groups.     ***/
  ReflectDVec(newtp->Charges, tp->Charges, tp->natom);
  for (i = 0; i < mp->nqvar-1; i++) {
    atmlist = mp->qvar[i].atoms;
    for (j = 0; j < mp->qvar[i].natm; j++) {
      newtp->Charges[atmlist[j]] = qval[i];
    }
  }
  atmlist = mp->qvar[i].atoms;
  for (j = 0; j < mp->qvar[i].natm; j++) {
    newtp->Charges[atmlist[j]] = 0.0;
  }
  qsum = -DSum(newtp->Charges, newtp->natom) / mp->qvar[i].natm;
  for (j = 0; j < mp->qvar[i].natm; j++) {
    newtp->Charges[atmlist[j]] = qsum;
  }

  /*** Free allocated memory ***/
  free(qval);

  /*** Increment counters ***/
  pivot = mp->nqvar-2;
  mp->qseed[pivot] += 1;
  while (mp->qseed[pivot] == mp->nqseed && pivot > 0) {
    mp->qseed[pivot] = 0;
    pivot--;
    mp->qseed[pivot] += 1;
  }
  if (pivot == 0 && mp->qseed[pivot] == mp->nqseed) {
    mp->qseed[0] = 0;
    return 1;
  }
  else {
    return 0;
  }

}

/***=======================================================================***/
/*** EvalParamFit: evaluate the degree to which a parameter set reproduces ***/
/***               the target energies.                                    ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:       fitting control data (contains history and number of      ***/
/***             fitting points)                                           ***/
/***   nrg:      the energies produced by this parameter set               ***/
/***   etrg:     the target energies                                       ***/
/***   ewt:      the weights of each conformation in the fit               ***/
/***=======================================================================***/
static double EvalParamFit(prmset *mp, double* nrg, double* etrg,
			   double* ewt)
{
  int i;
  double dE, wrmsd;

  wrmsd = 0.0;
  for (i = 0; i < mp->nconf; i++) {
    dE = nrg[i] - etrg[i];
    wrmsd += dE*dE*ewt[i];
  }

  return sqrt(wrmsd/mp->nconf);
}

/***=======================================================================***/
/*** ParamReport: report the best parameters and their fit to the target   ***/
/***              data.                                                    ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   tp:       the optimized topology                                    ***/
/***   etrg:     the target energies                                       ***/
/***   nrg:      the energies of all conformations according to the        ***/
/***             optimized topology                                        ***/
/***   mp:  fitting control data (contains the file name and the           ***/
/***             number of fitting points)                                 ***/
/***=======================================================================***/
static void ParamReport(double* etrg, double* nrg, prmset *mp, trajcon *tj)
{
  int i;
  FILE *outp;

  outp = FOpenSafe(tj->outbase, tj->OverwriteOutput);
  HorizontalRule(outp, 1);
  fprintf(outp, "   Target     Calculated\n   Energy       Energy  \n"
	  " ----------   ----------\n");
  for (i = 0; i < mp->nconf; i++) {
    fprintf(outp, " %10.6lf   %10.6lf\n", etrg[i], nrg[i]);
  }
  fprintf(outp, "\n");
  HorizontalRule(outp, 1);
  fclose(outp);
}

/***=======================================================================***/
/*** IncrementCharge: increment the charges belonging to one of the        ***/
/***                  variable groups in the molecule.                     ***/
/***=======================================================================***/
static void IncrementCharge(prmtop *tp, prmset *mp, int nq, double step)
{
  int i, j, nthisgrp, nothergrp;
  double qinc, qless;

  /*** Determine how many atoms are in this charge group ***/
  /*** versus how many atoms are in other charge groups. ***/
  nthisgrp = mp->qvar[nq].natm;
  nothergrp = 0;
  for (i = 0; i < mp->nqvar; i++) {
    if (i == nq) {
      nothergrp += mp->qvar[i].natm;
    }
  }
  qinc = step*(1.0e-5/nthisgrp);
  qless = -step*(1.0e-5/nothergrp);

  /*** Make the changes to the topology ***/
  for (i = 0; i < mp->qvar[nq].natm; i++) {
    tp->Charges[mp->qvar[nq].atoms[i]] += qinc;
  }
  for (i = 0; i < mp->nqvar; i++) {
    if (i == nq) {
      continue;
    }
    for (j = 0; j < mp->qvar[i].natm; j++) {
      tp->Charges[mp->qvar[i].atoms[j]] += qless;
    }
  }
}

/***=======================================================================***/
/*** OptimizeParam: optimize a given set of parameters with steepest       ***/
/***                descent minimization of the error relative to the      ***/
/***                target energies.                                       ***/
/***=======================================================================***/
static void OptimizeParam(prmtop *tp, imat *plist, dmat *crds, double* etrg,
			  double* ewt, prmset *mp)
{
  int i, j, k;
  double score, scorep, scoren, gmag, bestscore, newscore;
  double* nrg;
  double* bestnrg;
  double* grad;

  /*** Allocate space for current and best energy sets ***/
  nrg = (double*)malloc(mp->nconf*sizeof(double));
  bestnrg = (double*)malloc(mp->nconf*sizeof(double));
  grad = (double*)malloc(mp->nqvar*sizeof(double));
  for (i = 0; i < mp->maxiter; i++) {

    /*** Calculate the energies ***/
    ConfigNrg(tp, plist, crds, nrg);
    score = EvalParamFit(mp, nrg, etrg, ewt);

    /*** Compute the gradient in each variable charge group ***/
    gmag = 0.0;
    for (j = 0; j < mp->nqvar; j++) {

      /*** Get the +1 increment ***/
      IncrementCharge(tp, mp, j, 1.0);
      ConfigNrg(tp, plist, crds, nrg);
      scorep = EvalParamFit(mp, nrg, etrg, ewt);

      /*** Get the -1 increment ***/
      IncrementCharge(tp, mp, j, -2.0);
      ConfigNrg(tp, plist, crds, nrg);
      scoren = EvalParamFit(mp, nrg, etrg, ewt);

      /*** Return the topology to its unperturbed state ***/
      IncrementCharge(tp, mp, j, 1.0);
      grad[j] = scorep-scoren;
      gmag += grad[j]*grad[j];
    }

    /*** Normalize the gradient ***/
    gmag = 1.0/sqrt(gmag);
    for (j = 0; j < mp->nqvar; j++) {
      grad[j] *= gmag;
    }

    /*** Move down the gradient ***/
    printf("Energy %4d = %10.6lf\n", i, score);
    bestscore = score;
    for (j = 0; j < 100; j++) {
      for (k = 0; k < mp->nqvar; k++) {
	IncrementCharge(tp, mp, k, -grad[k]);
      }
      ConfigNrg(tp, plist, crds, nrg);
      newscore = EvalParamFit(mp, nrg, etrg, ewt);
      if (newscore < bestscore) {
	bestscore = newscore;
      }
      else {
	for (k = 0; k < mp->nqvar; k++) {
	  IncrementCharge(tp, mp, k, grad[k]);
	}
	j = 100;
      }
    }
  }

  /*** Free allocated memory ***/
  free(nrg);
  free(bestnrg);
  free(grad);
}

/***=======================================================================***/
/*** FitParams: the main function for optimizing parameters to fit a given ***/
/***            potential energy surface.                                  ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   mp:  fitting control data (contains the file name and the           ***/
/***             number of fitting points)                                 ***/
/***   tp:       the original topology, as read from a file                ***/
/***=======================================================================***/
void FitParams(prmset *mp, prmtop *tp, trajcon *tj)
{
  int i;
  double newscore, bestscore;
  double* etrg;
  double* nrg;
  double* bestnrg;
  double* ewt;
  imat plist;
  dmat crds;
  prmtop newtp;

  /*** Read the trajectory of conformations ***/
  crds = ReadSystemConf(mp, tp);

  /*** Read the list of energy targets ***/
  etrg = ReadTargetVal(mp, mp->sysnrg);
  ewt = ComputeNrgWeights(mp, etrg);

  /*** Make the pair list ***/
  plist = CreateNBList(tp);

  /*** Create new topologies to hold trial and best parameters ***/
  newtp = CopyTopology(tp);

  /*** Assign charge and Lennard-Jones groups ***/
  AssignChargeGroups(mp, tp);

  /*** Set charge and Lennard-Jones counters ***/
  mp->qseed = (int*)calloc(mp->nqvar, sizeof(int));
  mp->qval = (double*)malloc(mp->nqvar*sizeof(int));
  mp->bqval = (double*)malloc(mp->nqvar*sizeof(int));

  /*** Loop over all adjustable parameters and ***/
  /*** perform steepest-descent optimization   ***/
  nrg = (double*)malloc(mp->nconf*sizeof(double));
  bestnrg = (double*)malloc(mp->nconf*sizeof(double));
  for (i = 0; i < mp->nqseed; i++) {

    /*** Create the new Hamiltonian ***/
    SeedNewCharges(tp, &newtp, mp);

    /*** Optimize the new Hamiltonian ***/
    OptimizeParam(&newtp, &plist, &crds, etrg, ewt, mp);
    ConfigNrg(&newtp, &plist, &crds, nrg);
    newscore = EvalParamFit(mp, nrg, etrg, ewt);
    if (i == 0 || newscore < bestscore) {
      bestscore = newscore;
      ReflectDVec(bestnrg, nrg, mp->nconf);
    }
  }

  /*** Report the best parameters and the fit ***/
  ParamReport(etrg, bestnrg, mp, tj);

  /*** Free allocated memory ***/
  free(nrg);
  free(bestnrg);
  free(etrg);
  free(ewt);

  /*** Exit!  We are done. ***/
  exit(1);
}
