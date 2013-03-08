#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mdgxVector.h"
#include "Bonded.h"
#include "CrdManip.h"
#include "CellManip.h"
#include "Debug.h"

#include "TopologyDS.h"
#include "TrajectoryDS.h"

/***=======================================================================***/
/*** ElimInteraction: atoms X and Y may interact, but they should not      ***/
/***                  contribute to the total energy and forces because of ***/
/***                  exclusions.  This routine will test whether the      ***/
/***                  interaction of atoms X and Y has contributed to      ***/
/***                  energies, forces, or virials, then determine what    ***/
/***                  needs to be done to remove that contribution.        ***/
/***=======================================================================***/
static void ElimInteraction(coord *crd, nixpr* elist, int nelim, prmtop *tp,
			    int iordr, FrcTab *Cfrc, FrcTab *Hfrc,
			    Energy *sysUV)
{
  int i, xid, yid;
  double vdwscl, elecscl;
  double *xloc, *yloc;

  if (iordr == 4) {
    vdwscl = tp->lj14fac;
    elecscl = tp->elec14fac;
  }
  else {
    vdwscl = 1.0;
    elecscl = 1.0;
  }
  for (i = 0; i < nelim; i++) {
    xid = elist[i].atmX;
    yid = elist[i].atmY;
    xloc = &crd->scrloc[3*xid];
    yloc = &crd->scrloc[3*yid];
    if (sysUV->updateU == 2) {
      AttenuatePairVir(tp, Cfrc, Hfrc, elist[i].atmX, elist[i].atmY,
		       yloc[0]-xloc[0], yloc[1]-xloc[1], yloc[2]-xloc[2], 0.0,
		       &crd->scrfrc[3*xid], &crd->scrfrc[3*yid], sysUV,
		       elecscl, vdwscl);
    }
    else if (sysUV->updateU == 1) {
      AttenuatePairFrcNrg(tp, Cfrc, Hfrc, elist[i].atmX, elist[i].atmY,
		          yloc[0]-xloc[0], yloc[1]-xloc[1], yloc[2]-xloc[2],
			  0.0, &crd->scrfrc[3*xid], &crd->scrfrc[3*yid], sysUV,
			  elecscl, vdwscl);
    }
    else if (sysUV->updateU == 0) {
      AttenuatePairFrc(tp, Cfrc, Hfrc, elist[i].atmX, elist[i].atmY,
		       yloc[0]-xloc[0], yloc[1]-xloc[1], yloc[2]-xloc[2],
		       0.0, &crd->scrfrc[3*xid], &crd->scrfrc[3*yid], elecscl,
		       vdwscl);
    }
    else if (sysUV->updateU == -1) {
      AttenuatePairNrg(tp, Cfrc, Hfrc, elist[i].atmX, elist[i].atmY,
		       yloc[0]-xloc[0], yloc[1]-xloc[1], yloc[2]-xloc[2],
		       sysUV, elecscl, vdwscl);
    }
  }
}

/***=======================================================================***/
/*** CellBondedIntr: compute bonded interactions for a cell decomposition  ***/
/***                 of atom coordinates.  This routine can be called only ***/
/***                 after the "ShareTriple2" function in the CellManip    ***/
/***                 library has shared atoms between cells.  This routine ***/
/***                 searches for atoms within a region defined as the     ***/
/***                 same size as the home cell C, but translated by one   ***/
/***                 half the nonbonded cutoff in all directions.          ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   C:      the cell of interest                                        ***/
/***   CG:     the cell grid (used for dimensions)                         ***/
/***   tp:     the topology                                                ***/
/***   crd:    the coordinates (used for transformation matrices)          ***/
/***   Cfrc:   the "coarse" electrostatic force/energy spline table        ***/
/***   Hfrc:   the "fine" electrostatic force/energy spline table          ***/
/***   sysUV:  the system energy and virial (state information)            ***/
/***=======================================================================***/
void CellBondedIntr(cell *C, cellgrid *CG, coord *crd, prmtop *tp,
		    FrcTab *Cfrc, FrcTab *Hfrc, Energy *sysUV)
{
  int i, j, k, atmid, aid, bid, cid, did;
  double *atmcrd, *atmfrc;
  bondlist *tblc;
  bondcomm *tbcom;
  angllist *talc;
  anglcomm *tacom;
  dihelist *thlc;
  dihecomm *thcom;
  auxelim *eptmp;

  /*** Upload cell atom positons to the scratch array ***/
  UploadCellPosForce(C, crd);

  /*** Loop over all regions of this cell ***/
  for (i = 0; i < 8; i++) {

    /*** Loop over all atoms within the cell region ***/
    for (j = 0; j < C->nr[i]; j++) {
      if (IsCentralAtom(C->orig, &CG->celldim[4], C->map[i][j].loc,
			crd->U.data, CG->dbng, crd->isortho) == -1) {
	continue;
      }

      /*** If we are still here, then this atom may control ***/
      /*** bonded interactions that should be computed.     ***/
      atmid = C->map[i][j].id;
      atmcrd = &crd->scrloc[3*atmid];
      atmfrc = &crd->scrfrc[3*atmid];

      /*** Bonds ***/
      if (tp->BLC[atmid].nbond > 0) {
	tblc = &tp->BLC[atmid];
	for (k = 0; k < tblc->nbond; k++) {
	  tbcom = &tblc->BC[k];
	  bid = tbcom->b;
	  if (sysUV->updateU == 2) {
	    BondVir(atmcrd, &crd->scrloc[3*bid], atmfrc,
		    &crd->scrfrc[3*bid], tbcom, tp, Cfrc, Hfrc, sysUV);
	  }
	  else if (sysUV->updateU == 1) {
	    BondFrcNrg(atmcrd, &crd->scrloc[3*bid], atmfrc,
		       &crd->scrfrc[3*bid], tbcom, tp, Cfrc, Hfrc, sysUV);
	  }
	  else if (sysUV->updateU == -1) {
	    BondNrg(atmcrd, &crd->scrloc[3*bid], tbcom, tp, Cfrc, Hfrc,
		    sysUV);
	  }
	  else {
	    BondFrc(atmcrd, &crd->scrloc[3*bid], atmfrc, &crd->scrfrc[3*bid],
		    tbcom, tp, Cfrc, Hfrc);
	  }
	}
      }

      /*** Angles ***/
      if (tp->ALC[atmid].nangl > 0) {
	talc = &tp->ALC[atmid];
	for (k = 0; k < talc->nangl; k++) {
	  tacom = &talc->AC[k];
	  aid = tacom->a;
	  cid = tacom->c;
	  if (sysUV->updateU == 2) {
	    AnglVir(&crd->scrloc[3*aid], atmcrd, &crd->scrloc[3*cid],
		    &crd->scrfrc[3*aid], atmfrc, &crd->scrfrc[3*cid],
		    tacom, tp, Cfrc, Hfrc, sysUV);
	  }
	  else if (sysUV->updateU == 1) {
	    AnglFrcNrg(&crd->scrloc[3*aid], atmcrd,
		       &crd->scrloc[3*cid], &crd->scrfrc[3*aid],
		       atmfrc, &crd->scrfrc[3*cid], tacom, tp, Cfrc,
		       Hfrc, sysUV);
	  }
	  else if (sysUV->updateU == -1) {
	    AnglNrg(&crd->scrloc[3*aid], atmcrd, &crd->scrloc[3*cid],
		    tacom, tp, Cfrc, Hfrc, sysUV);
	  }
	  else {
	    AnglFrc(&crd->scrloc[3*aid], atmcrd, &crd->scrloc[3*cid],
		    &crd->scrfrc[3*aid], atmfrc, &crd->scrfrc[3*cid],
		    tacom, tp, Cfrc, Hfrc);
	  }
	}
      }

      /*** Dihedrals ***/
      if (tp->HLC[atmid].ndihe > 0) {
	thlc = &tp->HLC[atmid];
	for (k = 0; k < thlc->ndihe; k++) {
	  thcom = &thlc->HC[k];
	  aid = thcom->a;
	  cid = thcom->c;
	  did = thcom->d;
	  if (sysUV->updateU == 2) {
	    DiheVir(&crd->scrloc[3*aid], atmcrd, &crd->scrloc[3*cid],
		    &crd->scrloc[3*did], &crd->scrfrc[3*aid], atmfrc,
		    &crd->scrfrc[3*cid], &crd->scrfrc[3*did], thcom,
		    tp, Cfrc, Hfrc, sysUV);
	  }
	  else if (sysUV->updateU == 1) {
	    DiheFrcNrg(&crd->scrloc[3*aid], atmcrd,
		       &crd->scrloc[3*cid], &crd->scrloc[3*did],
		       &crd->scrfrc[3*aid], atmfrc,
		       &crd->scrfrc[3*cid], &crd->scrfrc[3*did], thcom,
		       tp, Cfrc, Hfrc, sysUV);
	  }
	  else if (sysUV->updateU == -1) {
	    DiheNrg(&crd->scrloc[3*aid], atmcrd, &crd->scrloc[3*cid],
		    &crd->scrloc[3*did], thcom, tp, Cfrc, Hfrc, sysUV);
	  }
	  else if (sysUV->updateU == 0) {
	    DiheFrc(&crd->scrloc[3*aid], atmcrd, &crd->scrloc[3*cid],
		    &crd->scrloc[3*did], &crd->scrfrc[3*aid], atmfrc,
		    &crd->scrfrc[3*cid], &crd->scrfrc[3*did], thcom,
		    tp, Cfrc, Hfrc);
	  }
	}
      }

      /*** If there have been no extra points added    ***/
      /*** to the topology and there are no additional ***/
      /*** exclusions to perform, we can just continue ***/
      if (tp->EPInserted == 0 && tp->ExclMarked == 0) {
	continue;
      }
      eptmp = &tp->ElimPair[atmid];
      ElimInteraction(crd, eptmp->list11, eptmp->n11, tp, 1, Cfrc, Hfrc,
		      sysUV);
      ElimInteraction(crd, eptmp->list12, eptmp->n12, tp, 2, Cfrc, Hfrc,
		      sysUV);
      ElimInteraction(crd, eptmp->list13, eptmp->n13, tp, 3, Cfrc, Hfrc,
		      sysUV);
      ElimInteraction(crd, eptmp->list14, eptmp->n14, tp, 4, Cfrc, Hfrc,
		      sysUV);
    }
  }

  DownloadCellForces(C, crd);
}

/*** BLOCKS 1 and 2: the forces are required ***/
#define NEEDFORCE 1

/*** BLOCK 1: energy (and virial) are not needed ***/
#define NEEDENERGY 0
#define NEEDVIRIAL 0
#include "BondContrib.c"
#undef NEEDVIRIAL
#undef NEEDENERGY

/*** BLOCK 2a: energy (but not virial) is needed, ***/
/***           in addition to forces              ***/
#define NEEDENERGY 1
#define NEEDVIRIAL 0
#include "BondContrib.c"
#undef NEEDVIRIAL

/*** BLOCK 2b: energy, virial, and forces are all needed ***/
#define NEEDVIRIAL 1
#include "BondContrib.c"
#undef NEEDVIRIAL
#undef NEEDENERGY

#undef NEEDFORCE
/*** END BLOCKS 1 and 2 ***/

/*** BLOCK 3: energy alone, but not forces or virials, are needed ***/
#define NEEDFORCE 0
#define NEEDENERGY 1
#define NEEDVIRIAL 0
#include "BondContrib.c"
#undef NEEDVIRIAL
#undef NEEDENERGY
#undef NEEDFORCE
/*** END BLOCK 3 ***/
