/***=======================================================================***/
/*** Cull2Point: cull atoms in a partitioned list according to their       ***/
/***             distance from either of two points.                       ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   btmp:     the list of atoms in this cell sector                     ***/
/***=======================================================================***/
#if VGTE == 1
static int Cull2PointVgtE(atomc *btmp, int *pordr, int natm, double* p1,
			  double* p2, double* pq1, double* pq2, int dim1,
			  int dim2, double Ecut2, double Vcut2)
#else
static int Cull2PointVeqE(atomc *btmp, int *pordr, int natm, double* p1,
			  double* p2, int dim1, int dim2, double Vcut2)
#endif
{
  int nact, i, poi;
  double r2, dv1, dv2;
  double *tloc;

  /*** Test and remove atoms which are not within bounds ***/
  nact = 0;
  for (i = 0; i < natm; i++) {
    poi = pordr[i];
    tloc = btmp[poi].loc;
#if VGTE == 1
    if (btmp[poi].lj >= 0) {
      dv1 = tloc[dim1] - p1[dim1];
      dv2 = tloc[dim2] - p1[dim2];
      r2 = dv1*dv1 + dv2*dv2;
      if (r2 < Vcut2) {
	pordr[nact] = poi;
	nact++;
      }
      else {
	dv1 = tloc[dim1] - p2[dim1];
	dv2 = tloc[dim2] - p2[dim2];
	r2 = dv1*dv1 + dv2*dv2;
	if (r2 < Vcut2) {
	  pordr[nact] = poi;
	  nact++;
	}
      }
    }
    else {
      dv1 = tloc[dim1] - pq1[dim1];
      dv2 = tloc[dim2] - pq1[dim2];
      r2 = dv1*dv1 + dv2*dv2;
      if (r2 < Ecut2) {
        pordr[nact] = poi;
        nact++;
      }
      else {
        dv1 = tloc[dim1] - pq2[dim1];
        dv2 = tloc[dim2] - pq2[dim2];
        r2 = dv1*dv1 + dv2*dv2;
        if (r2 < Ecut2) {
          pordr[nact] = poi;
          nact++;
        }
      }
    }
#else
    dv1 = tloc[dim1] - p1[dim1];
    dv2 = tloc[dim2] - p1[dim2];
    r2 = dv1*dv1 + dv2*dv2;
    if (r2 < Vcut2) {
      pordr[nact] = poi;
      nact++;
    }
    else {
      dv1 = tloc[dim1] - p2[dim1];
      dv2 = tloc[dim2] - p2[dim2];
      r2 = dv1*dv1 + dv2*dv2;
      if (r2 < Vcut2) {
        pordr[nact] = poi;
        nact++;
      }
    }
#endif
  }

  return nact;
}

/***=======================================================================***/
/*** Cull3Point: cull atoms in a partitioned list according to their       ***/
/***             distance from any of three points.                        ***/
/***=======================================================================***/
#if VGTE == 1
static int Cull3PointVgtE(atomc *btmp, int *pordr, int natm, double* p1,
			  double* p2, double* p3, double* pq1, double* pq2,
			  double* pq3, double Ecut2, double Vcut2)
#else
static int Cull3PointVeqE(atomc *btmp, int *pordr, int natm, double* p1,
			  double* p2, double* p3, double Vcut2)
#endif
{
  int nact, i, poi;
  double r2, dv1, dv2, dv3;
  double *tloc;

  /*** Test and remove atoms which are not within bounds ***/
  nact = 0;
  for (i = 0; i < natm; i++) {
    poi = pordr[i];
    tloc = btmp[poi].loc;
#if VGTE == 1
    if (btmp[poi].lj >= 0) {
      dv1 = tloc[0] - p1[0];
      dv2 = tloc[1] - p1[1];
      dv3 = tloc[2] - p1[2];
      r2 = dv1*dv1 + dv2*dv2 + dv3*dv3;
      if (r2 < Vcut2) {
	pordr[nact] = poi;
	nact++;
      }
      else {
	dv1 = tloc[0] - p2[0];
	dv2 = tloc[1] - p2[1];
	dv3 = tloc[2] - p2[2];
	r2 = dv1*dv1 + dv2*dv2 + dv3*dv3;
	if (r2 < Vcut2) {
	  pordr[nact] = poi;
	  nact++;
	}
	else {
	  dv1 = tloc[0] - p3[0];
	  dv2 = tloc[1] - p3[1];
	  dv3 = tloc[2] - p3[2];
	  r2 = dv1*dv1 + dv2*dv2 + dv3*dv3;
	  if (r2 < Vcut2) {
	    pordr[nact] = poi;
	    nact++;
	  }
	}
      }
    }
    else {
      dv1 = tloc[0] - pq1[0];
      dv2 = tloc[1] - pq1[1];
      dv3 = tloc[2] - pq1[2];
      r2 = dv1*dv1 + dv2*dv2 + dv3*dv3;
      if (r2 < Ecut2) {
        pordr[nact] = poi;
        nact++;
      }
      else {
        dv1 = tloc[0] - pq2[0];
        dv2 = tloc[1] - pq2[1];
        dv3 = tloc[2] - pq2[2];
        r2 = dv1*dv1 + dv2*dv2 + dv3*dv3;
        if (r2 < Ecut2) {
          pordr[nact] = poi;
          nact++;
        }
        else {
          dv1 = tloc[0] - pq3[0];
          dv2 = tloc[1] - pq3[1];
          dv3 = tloc[2] - pq3[2];
          r2 = dv1*dv1 + dv2*dv2 + dv3*dv3;
          if (r2 < Ecut2) {
            pordr[nact] = poi;
            nact++;
          }
        }
      }
    }
#else
    dv1 = tloc[0] - p1[0];
    dv2 = tloc[1] - p1[1];
    dv3 = tloc[2] - p1[2];
    r2 = dv1*dv1 + dv2*dv2 + dv3*dv3;
    if (r2 < Vcut2) {
      pordr[nact] = poi;
      nact++;
    }
    else {
      dv1 = tloc[0] - p2[0];
      dv2 = tloc[1] - p2[1];
      dv3 = tloc[2] - p2[2];
      r2 = dv1*dv1 + dv2*dv2 + dv3*dv3;
      if (r2 < Vcut2) {
	pordr[nact] = poi;
	nact++;
      }
      else {
	dv1 = tloc[0] - p3[0];
	dv2 = tloc[1] - p3[1];
	dv3 = tloc[2] - p3[2];
	r2 = dv1*dv1 + dv2*dv2 + dv3*dv3;
	if (r2 < Vcut2) {
	  pordr[nact] = poi;
	  nact++;
	}
      }
    }
#endif
  }

  return nact;
}

/***=======================================================================***/
/*** OrderOrtho: this function organizes the atoms in two sectors of a     ***/
/***             a cell structure.  Whether they share a face, edge, or    ***/
/***             just the center point is specified as an argument.        ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   C:      the cell whose regions are being ordered                    ***/
/***   cdim:   the dimensions of the cell's primary sector                 ***/
/***   vv:     modifier for the interaction type otype                     ***/
/***   otype:  the interaction type (0 = across face, 1 = across edge with ***/
/***           the primary sector, 2 = across edge with +Z sector, 3 =     ***/
/***           across edge with +Y sector, 4 = across cell center point)   ***/
/***   dcinp:  the direct space control information (for cutoffs)          ***/
/***=======================================================================***/
#if VGTE == 1
static void OrderOrthoVgtE(cell *C, double* cdim, int vv, int otype,
			   dircon *dcinp)
#else
static void OrderOrthoVeqE(cell *C, double* cdim, int vv, int otype,
			   dircon *dcinp)
#endif

{
  int i, j, natmA, natmB, maxatom, d1, d2;
  int *cnsr, *cordr, *cordr4, *itmp;
  double dx, dy, dz, dv, dr2, cmv1, cmv2;
  double cploc[3], mkAloc[3], mkBloc[3], mkCloc[3], cut[4];
#if VGTE == 1
  double qcut[4], mkAQloc[3], mkBQloc[3], mkCQloc[3];
#endif
  atomc *atmp, *btmp;

  /*** Determine sub-region cutoffs ***/
  cut[3] = dcinp->Mcut*dcinp->Mcut;
#if VGTE == 1
  qcut[3] = dcinp->Ecut*dcinp->Ecut;
#endif
  if (otype == 0) {
    cut[0] = 0.0625*cut[3];
    cut[1] = 0.2500*cut[3];
    cut[2] = 0.5625*cut[3];
#if VGTE == 1
    qcut[0] = 0.0625*qcut[3];
    qcut[1] = 0.2500*qcut[3];
    qcut[2] = 0.5625*qcut[3];
#endif
  }
  else {
    cut[0] = 0.25*cut[3];
    cut[1] = 0.50*cut[3];
    cut[2] = 0.75*cut[3];
#if VGTE == 1
    qcut[0] = 0.25*qcut[3];
    qcut[1] = 0.50*qcut[3];
    qcut[2] = 0.75*qcut[3];
#endif
  }

  /*** Location of the cell midpoint ***/
  for (i = 0; i < 3; i++) {
    cploc[i] = C->orig[i] + cdim[i];
  }

  /*** Unpack structures ***/
  cnsr = C->nsr;
  cordr = C->ordr.data;
  maxatom = C->maxatom;
  cordr4 = &cordr[4*maxatom];
  if (otype == 0) {

    /*** Interactions across a face ***/
    natmA = C->nr[0];
    atmp = C->map[0];
    if (vv == 0) {
      natmB = C->nr[1];
      btmp = C->map[1];
    }
    else if (vv == 1) {
      natmB = C->nr[2];
      btmp = C->map[2];
    }
    else if (vv == 2) {
      natmB = C->nr[4];
      btmp = C->map[4];
    }
  }
  else if (otype == 1) {

    /*** Interactions across an edge involving the primary sector ***/
    natmA = C->nr[0];
    atmp = C->map[0];
    if (vv == 2) {
      natmB = C->nr[3];
      btmp = C->map[3];
    }
    else if (vv == 1) {
      natmB = C->nr[5];
      btmp = C->map[5];
    }
    else if (vv == 0) {
      natmB = C->nr[6];
      btmp = C->map[6];
    }
  }
  else if (otype == 2) {

    /*** Interactions across an edge involving the +Z sector ***/
    natmA = C->nr[4];
    atmp = C->map[4];
    if (vv == 1) {
      natmB = C->nr[1];
      btmp = C->map[1];
    }
    else if (vv == 0) {
      natmB = C->nr[2];
      btmp = C->map[2];
    }
  }
  else if (otype == 3) {

    /*** Interactions across an edge involving the +Y sector ***/
    natmA = C->nr[2];
    atmp = C->map[2];
    if (vv == 2) {
      natmB = C->nr[1];
      btmp = C->map[1];
    }
  }
  else if (otype == 4) {

    /*** Interactions across the center point ***/
    if (vv == 0) {
      natmA = C->nr[0];
      atmp = C->map[0];
      natmB = C->nr[7];
      btmp = C->map[7];
    }
    else if (vv == 1) {
      natmA = C->nr[1];
      atmp = C->map[1];
      natmB = C->nr[6];
      btmp = C->map[6];
    }
    else if (vv == 2) {
      natmA = C->nr[2];
      atmp = C->map[2];
      natmB = C->nr[5];
      btmp = C->map[5];
    }
    else if (vv == 3) {
      natmA = C->nr[3];
      atmp = C->map[3];
      natmB = C->nr[4];
      btmp = C->map[4];
    }
  }

  /*** Prepare counters for the numbers of atoms in each sub-region ***/
  for (i = 0; i < 8; i++) {
    cnsr[i] = i*maxatom;
  }

  /*** First, do sub-region A ***/
  for (i = 0; i < natmA; i++) {

    /*** Compute the distance to the plane, edge, or point ***/
    if (otype == 0) {
      dr2 = cploc[vv] - atmp[i].loc[vv];
      dr2 *= dr2;
    }
    else {
      dx = cploc[0] - atmp[i].loc[0];
      dy = cploc[1] - atmp[i].loc[1];
      dz = cploc[2] - atmp[i].loc[2];
      dr2 = dx*dx + dy*dy + dz*dz;
      if (otype <= 3) {
        dv = cploc[vv] - atmp[i].loc[vv];
        dr2 -= dv*dv;
      }
    }

#if VGTE == 1
    /*** Eliminate atoms that have no van-der Waals properties ***/
    /*** if they are beyond the electrostatic cutoff.          ***/
    if (atmp[i].lj < 0) {
      if (dr2 < qcut[3]) {
	if (dr2 < qcut[1]) {
	  if (dr2 < qcut[0]) {
	    cordr[cnsr[0]] = i;
	    cnsr[0] += 1;
	  }
	  else {
	    cordr[cnsr[1]] = i;
	    cnsr[1] += 1;
	  }
	}
	else {
	  if (dr2 < qcut[2]) {
	    cordr[cnsr[2]] = i;
	    cnsr[2] += 1;
	  }
	  else {
	    cordr[cnsr[3]] = i;
	    cnsr[3] += 1;
	  }
	}
      }
      continue;
    }
#endif

    /*** Partition based on distance ***/
    if (dr2 < cut[3]) {
      if (dr2 < cut[1]) {
        if (dr2 < cut[0]) {
          cordr[cnsr[0]] = i;
          cnsr[0] += 1;
        }
        else {
          cordr[cnsr[1]] = i;
          cnsr[1] += 1;
        }
      }
      else {
        if (dr2 < cut[2]) {
          cordr[cnsr[2]] = i;
          cnsr[2] += 1;
        }
        else {
          cordr[cnsr[3]] = i;
          cnsr[3] += 1;
        }
      }
    }
  }

  /*** Next, do sub-region B ***/
  for (i = 0; i < natmB; i++) {

    /*** Compute the distance to the plane, edge, or point ***/
    if (otype == 0) {
      dr2 = cploc[vv] - btmp[i].loc[vv];
      dr2 *= dr2;
    }
    else {
      dx = cploc[0] - btmp[i].loc[0];
      dy = cploc[1] - btmp[i].loc[1];
      dz = cploc[2] - btmp[i].loc[2];
      dr2 = dx*dx + dy*dy + dz*dz;
      if (otype <= 3) {
        dv = cploc[vv] - btmp[i].loc[vv];
        dr2 -= dv*dv;
      }
    }

#if VGTE == 1
    /*** Eliminate atoms that have no van-der Waals properties ***/
    /*** if they are beyond the electrostatic cutoff.          ***/
    if (btmp[i].lj < 0) {
      if (dr2 < qcut[3]) {
	if (dr2 < qcut[1]) {
	  if (dr2 < qcut[0]) {
	    cordr[cnsr[4]] = i;
	    cnsr[4] += 1;
	  }
	  else {
	    cordr[cnsr[5]] = i;
	    cnsr[5] += 1;
	  }
	}
	else {
	  if (dr2 < qcut[2]) {
	    cordr[cnsr[6]] = i;
	    cnsr[6] += 1;
	  }
	  else {
	    cordr[cnsr[7]] = i;
	    cnsr[7] += 1;
	  }
	}
      }
      continue;
    }
#endif

    /*** Partition based on distance ***/
    if (dr2 < cut[3]) {
      if (dr2 < cut[1]) {
        if (dr2 < cut[0]) {
          cordr[cnsr[4]] = i;
          cnsr[4] += 1;
        }
        else {
          cordr[cnsr[5]] = i;
          cnsr[5] += 1;
        }
      }
      else {
        if (dr2 < cut[2]) {
          cordr[cnsr[6]] = i;
          cnsr[6] += 1;
        }
        else {
          cordr[cnsr[7]] = i;
          cnsr[7] += 1;
        }
      }
    }
  }

  /*** Adjust counters and compact the buffer regions ***/
  for (i = 1; i < 8; i++) {
    cnsr[i] -= i*maxatom;
  }
  for (i = 1; i < 4; i++) {
    cnsr[i] += cnsr[i-1];
    cnsr[i+4] += cnsr[i+3];
  }
  CompactOrder(cordr, cnsr[0], maxatom, cnsr[1]-cnsr[0]);
  CompactOrder(cordr, cnsr[1], 2*maxatom, cnsr[2]-cnsr[1]);
  CompactOrder(cordr, cnsr[2], 3*maxatom, cnsr[3]-cnsr[2]);
  CompactOrder(cordr4, cnsr[4], maxatom, cnsr[5]-cnsr[4]);
  CompactOrder(cordr4, cnsr[5], 2*maxatom, cnsr[6]-cnsr[5]);
  CompactOrder(cordr4, cnsr[6], 3*maxatom, cnsr[7]-cnsr[6]);

  /*** Copy the partitioned list of atoms in sub-region B ***/
  for (i = 1; i < 4; i++) {
    itmp = &cordr4[i*maxatom];
    for (j = 0; j < cnsr[7-i]; j++) {
      itmp[j] = cordr4[j];
    }
  }
  if (otype == 0) {
    return;
  }

  /*** If we're still here, we can cull a lot of atoms ***/
  for (i = 0; i < 3; i++) {
    cut[i] = sqrt(cut[i]);
#if VGTE == 1
    qcut[i] = sqrt(qcut[i]);
#endif
  }
  if (otype <= 3) {
    d1 = (vv == 0) ? 1 : 0;
    d2 = (vv == 2) ? 1 : 2;
    cmv1 = -1.0;
    cmv2 = (otype < 3) ? -1.0 : 1.0;
    for (i = 0; i < 3; i++) {
      mkAloc[d1] = cploc[d1] + cmv1*cut[i];
      mkAloc[d2] = cploc[d2];
      mkBloc[d1] = cploc[d1];
      mkBloc[d2] = cploc[d2] + cmv2*cut[i];
#if VGTE == 1
      mkAQloc[d1] = cploc[d1] + cmv1*qcut[i];
      mkAQloc[d2] = cploc[d2];
      mkBQloc[d1] = cploc[d1];
      mkBQloc[d2] = cploc[d2] + cmv2*qcut[i];
      cnsr[6-i] = Cull2PointVgtE(btmp, C->ordr.map[5+i], cnsr[6-i], mkAloc,
				 mkBloc, mkAQloc, mkBQloc, d1, d2, qcut[3],
				 cut[3]);
#else
      cnsr[6-i] = Cull2PointVeqE(btmp, C->ordr.map[5+i], cnsr[6-i], mkAloc,
				 mkBloc, d1, d2, cut[3]);
#endif
    }
  }

  if (otype == 4) {
    cmv1 = (vv == 0 || vv == 2) ? -1.0 : 1.0;
    cmv2 = (vv < 2) ? -1.0 : 1.0;
    for (i = 0; i < 3; i++) {
      mkAloc[0] = cploc[0] + cmv1*cut[i];
      mkAloc[1] = cploc[1];
      mkAloc[2] = cploc[2];
      mkBloc[0] = cploc[0];
      mkBloc[1] = cploc[1] + cmv2*cut[i];
      mkBloc[2] = cploc[2];
      mkCloc[0] = cploc[0];
      mkCloc[1] = cploc[1];
      mkCloc[2] = cploc[2] - cut[i];
#if VGTE == 1
      mkAQloc[0] = cploc[0] + cmv1*qcut[i];
      mkAQloc[1] = cploc[1];
      mkAQloc[2] = cploc[2];
      mkBQloc[0] = cploc[0];
      mkBQloc[1] = cploc[1] + cmv2*qcut[i];
      mkBQloc[2] = cploc[2];
      mkCQloc[0] = cploc[0];
      mkCQloc[1] = cploc[1];
      mkCQloc[2] = cploc[2] - qcut[i];
      cnsr[6-i] = Cull3PointVgtE(btmp, C->ordr.map[5+i], cnsr[6-i], mkAloc,
				 mkBloc, mkCloc, mkAQloc, mkBQloc, mkCQloc,
				 qcut[3], cut[3]);
#else
      cnsr[6-i] = Cull3PointVeqE(btmp, C->ordr.map[5+i], cnsr[6-i], mkAloc,
				 mkBloc, mkCloc, cut[3]);
#endif
    }
  }
}
