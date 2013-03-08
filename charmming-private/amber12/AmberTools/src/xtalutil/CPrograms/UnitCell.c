#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "pdbRead.h"
#include "crdmanip.h"
#include "matrix.h"
#include "myconstants.h"

struct symstruct {
  int nsym;
  dmat* tmat;
  double* tvec;
};
typedef struct symstruct symT;

symT LoadSymmetryData(char* source)
{
  int i, nm;
  double m1, m2, m3, t;
  char line[128], cword[128];
  FILE *inp;
  symT S;

  inp = fopen(source, "r");
  S.nsym = 0;
  while(fgets(line, 128, inp) != NULL) {
    if (strncmp(line, "REMARK", 6) == 0) {
      sscanf(line, "%s%s%s", cword, cword, cword);
      if (strncmp(cword, "SMTRY1", 6) == 0) {
	S.nsym += 1;
      }
    }
  }
  S.tmat = (dmat*)malloc(S.nsym*sizeof(dmat));
  for (i = 0; i < S.nsym; i++) {
    S.tmat[i] = CreateDmat(3, 3);
  }
  S.tvec = (double*)malloc(3*S.nsym*sizeof(double));
  rewind(inp);
  while(fgets(line, 128, inp) != NULL) {
    if (strncmp(line, "REMARK", 6) == 0) {
      sscanf(line, "%s%s%s", cword, cword, cword);
      if (strncmp(cword, "SMTRY", 5) == 0) {
	sscanf(line, "%s%s%s%d%lf%lf%lf%lf\n", cword, cword, cword, &nm, &m1,
	       &m2, &m3, &t);
	cword[6] = '\0';
	i = atoi(&cword[5])-1;
	nm -= 1;
	S.tmat[nm].data[3*i] = m1;
	S.tmat[nm].data[3*i+1] = m2;
	S.tmat[nm].data[3*i+2] = m3;
	S.tvec[3*nm+i] = t;
      }
    }
  }
  fclose(inp);

  return S;
}
	      
int main(int argc, char *argv[])
{
  int h, i, m, n, tm, rskip;
  double x, y, z;
  double *dtmp;
  char tag[64];
  pdb p, t;
  symT S;

  if (argc == 1) {
    printf("\nUnitCell >> A program for recreating a crystallographic unit "
	   "cell from a PDB\n""UnitCell >> structure.\n\n"
	   "Options:\n"
           "  -p   : the structure to reassemble (PDB format)\n"
           "  -o   : the output structure (PDB format)\n");
    exit(1);
  }

  /*** Input ***/
  p.source[0] = '\0';
  t.source[0] = '\0';
  for (i = 0; i < argc-2; i += 2) {
    strcpy(tag, *++argv);
    if (strcmp(tag, "-p") == 0) {
      strcpy(p.source, *++argv);
    }
    else if (strcmp(tag, "-o") == 0) {
      strcpy(t.source, *++argv);
    }
    else {
      printf("UnitCell >> Error.  Tag %s not recognized.\n", tag);
      exit(1);
    }
  }

  /*** Check input ***/
  if (p.source[0] == '\0') {
    printf("UnitCell >> Error.  Original PDB file not specified!\n");
    exit(1);
  }
  if (t.source[0] == '\0') {
    printf("UnitCell >> Error.  Output PDB file not specified!\n");
    exit(1);
  }

  /*** Get the original PDB ***/
  GetPDB(&p, 0, 0);

  /*** Scan PDB file for symmetry information ***/
  S = LoadSymmetryData(p.source);

  /*** Now, make the new PDB ***/
  rskip = p.n_res;
  t.n_atoms = S.nsym*p.n_atoms;
  t.n_res = S.nsym*p.n_res;
  t.n_ter = S.nsym*(p.n_ter+1);
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

  /*** Construct symmetry-related copies ***/
  h = 0;
  for (i = 0; i < S.nsym; i++) {
    for (m = 0; m < p.n_atoms; m++) {
      tm = m + h*p.n_atoms;
      t.res_nums[tm] = p.res_nums[m] + h*rskip;
      t.atom_nums[tm] = p.atom_nums[m] + h*p.n_atoms;
      t.chain[tm] = p.chain[m];
      for (n = 0; n < 4; n++) {
	t.res_names[4*tm+n] = p.res_names[4*m+n];
	t.atom_names[4*tm+n] = p.atom_names[4*m+n];
      }
      dtmp = S.tmat[i].data;
      x = p.crds[3*m];
      y = p.crds[3*m+1];
      z = p.crds[3*m+2];
      t.crds[3*tm] = dtmp[0]*x + dtmp[1]*y + dtmp[2]*z + S.tvec[3*i];
      t.crds[3*tm+1] = dtmp[3]*x + dtmp[4]*y + dtmp[5]*z + S.tvec[3*i+1];
      t.crds[3*tm+2] = dtmp[6]*x + dtmp[7]*y + dtmp[8]*z + S.tvec[3*i+2];
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
  t.res_lims[t.n_res] = t.n_atoms;

  /*** Print PDB ***/
  ModPdbRA(&t);
  PutPDB(&t, t.source, "STANDARD", "HEADER  ", 0);

  return 0;
}
