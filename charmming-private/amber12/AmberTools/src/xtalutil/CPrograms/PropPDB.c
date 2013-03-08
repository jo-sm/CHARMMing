#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "pdbRead.h"
#include "crdmanip.h"
#include "matrix.h"
#include "myconstants.h"

int main(int argc, char *argv[])
{
  int h, i, j, k, m, n, tm, nrep, rskip;
  int irep[3];
  double boxd[6];
  char tag[64];
  dmat U, invU;
  pdb p, t;

  if (argc == 1) {
    printf("\nPropPDB >> A program for propagating a PDB structure.\n\n"
	   "Options:\n"
           "  -p   : the structure to reassemble (PDB format)\n"
           "  -o   : the output structure (PDB format)\n"
	   "  -X   : the length of the first unit cell vector\n"
	   "  -Y   : the length of the second unit cell vector\n"
	   "  -Z   : the length of the third unit cell vector\n"
	   "  -a   : the box alpha angle (degrees)\n"
	   "  -b   : the box beta angle (degrees)\n"
	   "  -g   : the box gamma angle (degrees)\n" 
	   "  -ix  : number of replicas along _X_ vector\n"
	   "  -iy  : number of replicas along _Y_ vector\n"
	   "  -iz  : number of replicas along _Z_ vector\n\n");
    exit(1);
  }

  /*** Input ***/
  p.source[0] = '\0';
  t.source[0] = '\0';
  boxd[0] = -1.0;
  boxd[1] = -1.0;
  boxd[2] = -1.0;
  boxd[3] = 90.0;
  boxd[4] = 90.0;
  boxd[5] = 90.0;
  irep[0] = 1;
  irep[1] = 1;
  irep[2] = 1;
  for (i = 0; i < argc-2; i += 2) {
    strcpy(tag, *++argv);
    if (strcmp(tag, "-p") == 0) {
      strcpy(p.source, *++argv);
    }
    else if (strcmp(tag, "-o") == 0) {
      strcpy(t.source, *++argv);
    }
    else if (strcmp(tag, "-X") == 0) {
      boxd[0] = atof(*++argv);
    }
    else if (strcmp(tag, "-Y") == 0) {
      boxd[1] = atof(*++argv);
    }
    else if (strcmp(tag, "-Z") == 0) {
      boxd[2] = atof(*++argv);
    }
    else if (strcmp(tag, "-a") == 0) {
      boxd[3] = atof(*++argv);
    }
    else if (strcmp(tag, "-b") == 0) {
      boxd[4] = atof(*++argv);
    }
    else if (strcmp(tag, "-g") == 0) {
      boxd[5] = atof(*++argv);
    }
    else if (strcmp(tag, "-ix") == 0) {
      irep[0] = atoi(*++argv);
    }
    else if (strcmp(tag, "-iy") == 0) {
      irep[1] = atoi(*++argv);
    }
    else if (strcmp(tag, "-iz") == 0) {
      irep[2] = atoi(*++argv);
    }
    else {
      printf("PropPDB >> Error.  Unrecognized tag %s.\n", tag);
      exit(1);
    }
  }

  /*** Check ***/
  if (p.source[0] == '\0') {
    printf("PropPDB >> Error.  Original PDB file not specified!\n");
    exit(1);
  }
  if (t.source[0] == '\0') {
    printf("PropPDB >> Error.  Output PDB file not specified!\n");
    exit(1);
  }
  if (boxd[0] < 0.0 || boxd[1] < 0.0 || boxd[2] < 0.0) {
    printf("PropPDB >> Error.  Box dimensions invalid!\n");
    exit(1);
  }

  /*** Transform the box ***/
  for (i = 0; i < 3; i++) {
    boxd[3+i] *= PI/180.0;
  }

  /*** Get the first PDB ***/
  GetPDB(&p, 1, 0);

  /*** Now, make the new PDB ***/
  nrep = irep[0]*irep[1]*irep[2];
  rskip = p.n_res;
  t.n_atoms = nrep*p.n_atoms;
  t.n_res = nrep*p.n_res;
  t.n_ter = nrep*(p.n_ter+1);
  t.aug = 0;
  t.pqr = 0;
  t.res_nums = (int*)malloc(t.n_atoms*sizeof(int));
  t.atom_nums = (int*)malloc(t.n_atoms*sizeof(int));
  t.res_lims = (int*)malloc((t.n_res+1)*sizeof(int));
  t.termini = (int*)malloc(t.n_ter*sizeof(int));
  t.chain = (char*)malloc(t.n_atoms*sizeof(char));
  t.atom_names = (char*)malloc(4*t.n_atoms*sizeof(char));
  t.res_names = (char*)malloc(4*t.n_atoms*sizeof(char));
  t.crds = (double*)malloc(3*t.n_atoms*sizeof(double));

  /*** Go into box space ***/
  U = CreateDmat(3, 3);
  invU = CreateDmat(3, 3);
  CmpXfrm(boxd, U, invU);
  RotateCrd(p.crds, p.n_atoms, U);

  /*** Propagate ***/
  for (m = 0; m < p.n_atoms; m++) {
    t.res_nums[m] = p.res_nums[m]; 
    t.atom_nums[m] = p.atom_nums[m]; 
    t.chain[m] = p.chain[m];
    for (n = 0; n < 4; n++) {
      t.res_names[4*m+n] = p.res_names[4*m+n];
      t.atom_names[4*m+n] = p.atom_names[4*m+n];
    }
    for (n = 0; n < 3; n++) {
      t.crds[3*m+n] = p.crds[3*m+n];
    }
  }
  for (m = 0; m < p.n_res; m++) {
    t.res_lims[m] = p.res_lims[m];
  }
  for (m = 0; m < p.n_ter; m++) {
    t.termini[m] = p.termini[m];
  }
  t.termini[m] = p.n_atoms;

  h = 0;
  for (i = 0; i < irep[0]; i++) {
    for (j = 0; j < irep[1]; j++) {
      for (k = 0; k < irep[2]; k++) {
	for (m = 0; m < p.n_atoms; m++) {
	  tm = m + h*p.n_atoms;
	  t.res_nums[tm] = p.res_nums[m] + h*rskip;
	  t.atom_nums[tm] = p.atom_nums[m] + h*p.n_atoms;
	  t.chain[tm] = p.chain[m];
	  for (n = 0; n < 4; n++) {
	    t.res_names[4*tm+n] = p.res_names[4*m+n];
	    t.atom_names[4*tm+n] = p.atom_names[4*m+n];
	  }
	  t.crds[3*tm] = p.crds[3*m] + i;
	  t.crds[3*tm+1] = p.crds[3*m+1] + j;
	  t.crds[3*tm+2] = p.crds[3*m+2] + k;
	}
	for (m = 0; m < p.n_res; m++) {
	  tm = m + h*p.n_res;
	  t.res_lims[tm] = p.res_lims[m] + h*p.n_atoms;
	}
	for (m = 0; m < p.n_ter; m++) {
	  tm = m + h*(p.n_ter+1);
	  t.termini[tm] = p.termini[m] + h*p.n_atoms;
	}
	t.termini[m + h*(p.n_ter+1)] = (h+1)*p.n_atoms;
	h++;
      }
    }
  }
  t.res_lims[t.n_res] = t.n_atoms;

  /*** Back to real space ***/
  RotateCrd(t.crds, t.n_atoms, invU);

  /*** Print PDB ***/
  ModPdbRA(&t);
  PutPDB(&t, t.source, "STANDARD", "HEADER  ", 1);

  return 0;
}
