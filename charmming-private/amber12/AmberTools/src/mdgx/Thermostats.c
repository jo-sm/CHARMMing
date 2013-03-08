#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "Random.h"
#include "mdgxVector.h"
#include "Thermostats.h"
#include "CellManip.h"

#include "TrajectoryDS.h"
#include "TopologyDS.h"
#include "CrdManipDS.h"

#define NEED_TI 0
#include "ThermostatsBranch.c"
#undef NEED_TI

#define NEED_TI 1
#include "ThermostatsBranch.c"
#undef NEED_TI

/***=======================================================================***/
/*** CellNHThermostat: a compact function for a nasty series of            ***/
/***                   expressions.  This implements the repeatedly used   ***/
/***                   thermostat in the Nose-Hoover thermo- (baro-)stat   ***/
/***                   above.  Note that if there is no barsotat (barostat ***/
/***                   is 0), then this works on a 1/4 timestep interval.  ***/
/***                   If there is, it works on a 1/8 timestep interval.   ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   CG:       the cell grid                                             ***/
/***   crd:      the coordinates (cells reference coordinates for velocity ***/
/***             information)                                              ***/
/***   tp:       the topology                                              ***/
/***   tj:       trajectory control information                            ***/
/***   sysUV:    the system energy and virial                              ***/
/***   ETAt:     scaling factor relating to the barostat piston mass       ***/
/***   CHIt:     (modified and returned) exponential argument for velocity ***/
/***               rescaling, similar to chi in the Berendsen thermostat   ***/
/***   CHItp1e:  CHIt advanced 1/8 timestep                                ***/
/***   CHItp1q:  CHIt advanced 1/4 timestep                                ***/
/***   barostat: flag to indicate whether a Nose-Hoover barostat is active ***/
/***=======================================================================***/
void NHThermostat(cellgrid *CG, coord *crd, prmtop *tp, trajcon *tj,
		  Energy *sysUV, double ETAt, double *CHIt, double *CHItp1e,
		  double *CHItp1q, int barostat)
{
  int i, j, g3con;
  double dt, Ttarget, pmass, qmass, invqmass, sigma, vfac;
  double *vtmp;
  cell *C;

  /*** Unpack ***/
  dt = (barostat == 1) ? 0.125*tj->dt : 0.25*tj->dt;
  Ttarget = (barostat == 1) ? GASCNST*tj->Ttarget : 0.0;
  pmass = (barostat == 1) ? ETAt*ETAt*tj->npth.pmass : 0.0;
  qmass = tj->npth.qmass;
  sigma = 2.0*tj->npth.sigma;
  invqmass = 1.0/qmass;

  /*** Compute the kinetic energy ***/
  sysUV->kine = KineticEnergy(CG, crd, tp);
  *CHItp1e = *CHIt + dt*(2.0*sysUV->kine + pmass - sigma - Ttarget)*invqmass;
  const int ncell = CG->ng[0]*CG->ng[1]*CG->ng[2];
  vtmp = crd->vel;
  vfac = exp(-2.0*dt*(*CHItp1e));
  for (i = 0; i < ncell; i++) {
    C = &CG->data[i];
    for (j = 0; j < C->nr[0]; j++) {
      g3con = 3*C->data[j].id;
      vtmp[g3con] *= vfac;
      vtmp[g3con+1] *= vfac;
      vtmp[g3con+2] *= vfac;
    }
  }
  sysUV->kine = KineticEnergy(CG, crd, tp);
  *CHItp1q = *CHItp1e +
    dt*(2.0*sysUV->kine + pmass - sigma - Ttarget)*invqmass;
}

#if 0
/***=======================================================================***/
/*** NPTHBarostat: a compact function for another nasty series of          ***/
/***               expressions in the coupled Nose-Hoover thermo-barostat  ***/
/***               encoded above.                                          ***/
/***=======================================================================***/
void NPTHBarostat(coord *crd, trajcon *tj, double CHIt, double *ETAt,
                  double *ETAtp1q, double *ETAtp1h)
{
  double CHIfac, dt, Vt, Ptarget, invpmass, invqmass;

  /*** Unpack ***/
  dt = rp->dt;
  Ptarget = rp->Ptarget;
  Vt = rp->stv.V;
  invpmass = 1.0/tj->npth.pmass;
  invqmass = 1.0/tj->npth.qmass;

  CHIfac = exp(-0.125*CHIt*dt);
  *ETAt *= CHIfac;

  rp->stv.P = InstantSystemPressure(rp);
  *ETAtp1q = *ETAt + 0.25*dt*(3.0*(rp->stv.P - Ptarget)*Vt)*invpmass;
  *ETAtp1q *= CHIfac;

  DVecMult(ma->vel, 3*ma->natm, exp(-0.5*dt*(*ETAtp1q)));
  rp->nrg.ke = KineticEnergyComp(ma);
  *ETAtp1q *= CHIfac;

  rp->stv.P = InstantSystemPressure(rp);
  *ETAtp1h = *ETAtp1q + 0.25*dt*(3.0*(rp->stv.P - Ptarget)*Vt)*invpmass;
}
#endif

/***=======================================================================***/
/*** PrepThermoBarostat: prepare constants for a thermostat and barostat,  ***/
/***                     if needed.                                        ***/
/***=======================================================================***/
void PrepThermoBarostat(prmtop *tp, trajcon *tj)
{
  tj->npth.sigma = 0.5*tp->ndf*GASCNST*tj->Ttarget;
  tj->npth.qmass = 2.0*tj->npth.sigma*tj->npth.TauT*tj->npth.TauT;
  tj->npth.pmass = (tp->ndf + 3)*GASCNST*tj->Ttarget*pow(tj->npth.TauP, 2.0);
  tj->npth.chi = 0.0;
  tj->npth.eta = 0.0;
}
