#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "pdbRead.h"

/***=======================================================================***/
/*** GetPDB: function for taking in the information within a standard-     ***/
/***         format PDB file.  The information is returend within a pdb    ***/
/***         type structure.                                               ***/
/***=======================================================================***/
void GetPDB(pdb *thispdb, int verbosity, int pqr_format)
{
  int i, j, t, curr_res, rnumrescue, anumrescue;
  char line[MAXLINE];
  FILE *pdbsource;

  /*** Scan the PDB file for atoms ***/
  if ((pdbsource = fopen(thispdb->source, "r")) == NULL) {
    printf("GetPDB >> Error.  PDB File %s not found!\n",
	   thispdb->source);
    exit(1);
  }
  thispdb->n_atoms = 0;
  thispdb->n_ter = 0;
  while(fgets(line, MAXLINE, pdbsource) != NULL) {
    if ((line[0] == 'A' && strncmp(line, "ATOM", 4) == 0) ||
	(line[0] == 'H' && strncmp(line, "HETATM", 6) == 0)) {
      thispdb->n_atoms += 1;
    }
    else if (strncmp(line, "TER", 3) == 0) {
      thispdb->n_ter += 1;
    }
  }
  if (verbosity == 1) {
    printf("GetPDB >> Found %d atoms in %s.\n", thispdb->n_atoms,
	   thispdb->source);
  }
  rewind(pdbsource);

  /*** Allocate memory for the molecule ***/
  thispdb->res_nums = (int*)malloc(thispdb->n_atoms*sizeof(int));
  thispdb->atom_nums = (int*)malloc(thispdb->n_atoms*sizeof(int));
  thispdb->chain = (char*)malloc(thispdb->n_atoms*sizeof(char));
  thispdb->atom_names = (char*)malloc(4*thispdb->n_atoms*sizeof(char));
  thispdb->res_names = (char*)malloc(4*thispdb->n_atoms*sizeof(char));
  thispdb->crds = (double*)calloc(3*thispdb->n_atoms, sizeof(double));
  thispdb->termini = (int*)malloc(thispdb->n_ter*sizeof(int));
  thispdb->pqr = pqr_format;
  thispdb->aug = 0;
  if (thispdb->pqr == 1) {
    thispdb->charges = (double*)malloc(thispdb->n_atoms*sizeof(double));
    thispdb->radii = (double*)malloc(thispdb->n_atoms*sizeof(double));
  }

  /*** Read in the PDB information ***/
  i = 0;
  t = 0;
  rnumrescue = 0;
  anumrescue = 0;
  while(fgets(line, MAXLINE, pdbsource) != NULL) {
    if ((line[0] == 'A' && strncmp(line, "ATOM", 4) == 0) ||
	(line[0] == 'H' && strncmp(line, "HETATM", 6) == 0)) {
      if (line[6] != '*') {
	sscanf(&line[4], "%d", &thispdb->atom_nums[i]);
	anumrescue = thispdb->atom_nums[i];
      }
      else {
	thispdb->atom_nums[i] = anumrescue;
	anumrescue++;
      }
      for (j = 0; j < 4; j++) {
	thispdb->atom_names[i*4+j] = line[j+12];
	thispdb->res_names[i*4+j] = line[j+17];
      }
      if (line[22] != '*') {
        if (thispdb->pqr == 0) {
          sscanf(&line[21], "%c%d%lf%lf%lf", &thispdb->chain[i],
	         &thispdb->res_nums[i], &thispdb->crds[i*3],
	         &thispdb->crds[i*3+1], &thispdb->crds[i*3+2]);
        }
        else {
          sscanf(&line[21], "%c%d%lf%lf%lf%lf%lf", &thispdb->chain[i],
	         &thispdb->res_nums[i], &thispdb->crds[i*3],
	         &thispdb->crds[i*3+1], &thispdb->crds[i*3+2],
	         &thispdb->charges[i], &thispdb->radii[i]);
	}
	rnumrescue = thispdb->res_nums[i]+1;
      }
      else {
	thispdb->chain[i] = line[21];
        if (thispdb->pqr == 0) {
          sscanf(&line[30], "%lf%lf%lf", &thispdb->crds[i*3],
	         &thispdb->crds[i*3+1], &thispdb->crds[i*3+2]);
        }
        else {
          sscanf(&line[30], "%lf%lf%lf%lf%lf", &thispdb->crds[i*3],
	         &thispdb->crds[i*3+1], &thispdb->crds[i*3+2],
	         &thispdb->charges[i], &thispdb->radii[i]);
	}
	thispdb->res_nums[i] = rnumrescue;
      }

      /*** Filter the character inputs ***/
      touchup(&thispdb->atom_names[i*4], 4);
      touchup(&thispdb->res_names[i*4], 4); alphanumeric(&thispdb->chain[i]);
      i++;
    }
    if (strncmp(line, "TER", 3) == 0) {
      thispdb->termini[t] = i;
      t++;
      rnumrescue++;
    }
  }
  fclose(pdbsource);

  /*** Determine the residue limits and normalize residue indices ***/
  thispdb->n_res = 1;
  curr_res = thispdb->res_nums[0];
  j = 0;
  for (i = 0; i < thispdb->n_atoms; i++) {
    if (thispdb->res_nums[i] != curr_res) {
      thispdb->n_res += 1;
      curr_res = thispdb->res_nums[i];
      j++;
    }
    thispdb->res_nums[i] = j;
  }
  thispdb->res_lims = (int*)malloc((thispdb->n_res+1)*sizeof(int));
  curr_res = thispdb->res_nums[0];
  thispdb->res_lims[0] = 0;
  j = 1;
  for (i = 0; i < thispdb->n_atoms; i++) {
    if (thispdb->res_nums[i] != curr_res) {
      curr_res = thispdb->res_nums[i];
      thispdb->res_lims[j] = i;
      j++;
    }
  }
  thispdb->res_lims[thispdb->n_res] = thispdb->n_atoms;

  /*** Make sure TER cards fall in between residues ***/
  for (i = 0; i < thispdb->n_ter; i++) {
    t = 0;
    for (j = i; j < thispdb->n_res+1; j++) {
      if (thispdb->res_lims[j] == thispdb->termini[i]) {
	t = 1;
	break;
      }
    }
    if (t == 0 && verbosity == 1) {
      printf("GetPDB >> Warning.  TER card after atom %d does not match "
	     "residue limits!\n", thispdb->termini[i]);
    }
  }

  /*** No TER card after the last residue ***/
  if (thispdb->n_ter > 0 &&
      thispdb->termini[thispdb->n_ter-1] == thispdb->n_atoms) {
    thispdb->n_ter -= 1;
  }

  if (verbosity == 1) {
    printf("GetPDB >> Bracketed %d residues.\n\n", thispdb->n_res);
  }
}

/***=======================================================================***/
/*** PutPDB: function for writing PDB files in a specified format, an      ***/
/***         improved version of the original and a companion to GetPDB.   ***/
/***=======================================================================***/
void PutPDB(pdb *thispdb, char* filename, char* pdb_form,
	    char* custom_header, int verbosity)
{
  int i, j, k, minj;
  double* dtmp;
  char* ctmp;
  char* ctm2p;
  char c_holding;
  char line[MAXLINE];
  FILE *pdbtxt;

  pdbtxt = fopen(filename, "w");
  fprintf(pdbtxt, "HEADER  Generated by GRAPPLE (D. Cerutti & T. Jain, "
	  "McCammon Lab, UCSD)\n%s\n", custom_header);
  minj = 0;
  for (i = 0; i < thispdb->n_atoms; i++) {
    ctmp = &thispdb->atom_names[i*4];
    ctm2p = &thispdb->res_names[i*4];
    dtmp = &thispdb->crds[i*3];

    /*** The basic information ***/
    if (thispdb->pqr == 0) {
      sprintf(line, "ATOM %6d %.4s %.4s%c%4d    %8.3lf%8.3lf%8.3lf\n",
	      thispdb->atom_nums[i], ctmp, ctm2p, thispdb->chain[i],
	      thispdb->res_nums[i], dtmp[0], dtmp[1], dtmp[2]);
    }
    else {
      sprintf(line, "ATOM %6d %.4s %.4s %4d    %8.3lf%8.3lf%8.3lf "
	      "%9.4lf %9.4lf\n", thispdb->atom_nums[i], ctmp, ctm2p,
	      thispdb->res_nums[i], dtmp[0], dtmp[1], dtmp[2],
	      thispdb->charges[i], thispdb->radii[i]);
    }

    /*** Print the line with modifications ***/
    if (strncmp(pdb_form, "STANDARD", 7) == 0) {

      /*** Standard format ***/
      realign_name(&line[12]);
      fprintf(pdbtxt, "%s", line);
    }
    else if (strncmp(pdb_form, "UHBD", 4) == 0) {

      /*** The atom names must begin with letters ***/
      ctmp = &line[12];
      realign_name(ctmp);
      if (ctmp[0] >= '0'  && ctmp[0] <= '9') {
	j = 3;
	while (j >= 0 && ctmp[j] == ' ') {
	  j--;
	}
	c_holding = ctmp[0];
        for (k = 0; k < j; k++) {
	  ctmp[k] = ctmp[k+1];
        }
	ctmp[j] = c_holding;
      }
      realign_name(&line[12]);
      fprintf(pdbtxt, "%s", line);
    }

    /*** Print TER card ***/
    for (j = 0; j < thispdb->n_ter; j++) {
      if (thispdb->termini[j] == i+1) {
	minj = j;
	fprintf(pdbtxt, "TER\n");
	break;
      }
      if (thispdb->termini[j] > i+1) {
	break;
      }
    }
  }
  fprintf(pdbtxt, "END\n");
  fclose(pdbtxt);
  if (verbosity == 1) {
    printf("PutPDB >> Molecule printed to file %s in %s format\n\n",
	   filename, pdb_form);
  }
}

/***=======================================================================***/
/*** FreePDB: function for freeing a pdb structure.  Easy as pie!          ***/
/***=======================================================================***/
void FreePDB(pdb *thispdb)
{
  free(thispdb->res_nums);
  free(thispdb->atom_nums);
  free(thispdb->res_lims);
  free(thispdb->chain);
  free(thispdb->atom_names);
  free(thispdb->res_names);
  free(thispdb->crds);
  free(thispdb->termini);
  if (thispdb->pqr == 1) {
    free(thispdb->charges);
    free(thispdb->radii);
  }
  if (thispdb->aug == 1) {
    free(thispdb->z);
    free(thispdb->mass);
  }
  thispdb->pqr = 0;
  thispdb->aug = 0;
}

/***=======================================================================***/
/*** SSCAN_ADV: scans a string from a larger string, then advances to the  ***/
/***            last character scanned in the larger line.  Ignores any    ***/
/***            whitespace.                                                ***/
/***=======================================================================***/
char* sscan_adv(char *ctmp, char* cword)
{
  while (ctmp[0] == ' ') {
    ctmp = &ctmp[1];
  }
  sscanf(ctmp, "%s", cword);
  ctmp = &ctmp[strlen(cword)];

  return ctmp;
}

/***=======================================================================***/
/*** TOUCHUP: combines the functions append_whitespace, remove_whitespace, ***/
/***          and alphanumeric into one convenient package.                ***/
/***=======================================================================***/
void touchup(char* a, int asize)
{
  int i;

  append_whitespace(a, asize);
  remove_whitespace(a, asize);
  for (i = 0; i < asize; i++) {
    if (a[i] == '+' || a[i] == '-') {
      continue;
    }
    alphanumeric(&a[i]);
  }
}

/***=======================================================================***/
/*** REMOVE_WHITESPACE: function for removing white space at the beginning ***/
/***                    of a string.                                       ***/
/***=======================================================================***/
void remove_whitespace(char* a, int asize)
{
  int i, j, num_rem;

  num_rem = 0;
  while (a[num_rem] == ' ' && num_rem < asize) {
    num_rem++;
  }
  j = 0;
  for (i = num_rem; i < asize; i++) {
    a[j] = a[i];
    j++;
  }
  for (i = asize-num_rem; i < asize; i++) {
    a[i] = ' ';
  }
}

/***=======================================================================***/
/*** APPEND_WHITESPACE: function for adding whitespace to the end of a     ***/
/***                    string, staring at the (first) instance of '\0'    ***/
/***=======================================================================***/
void append_whitespace(char* a, int asize)
{
  int i, whiteout;
  
  whiteout = 0;
  for (i = 0; i < asize; i++) {
    if (a[i] == '\0') {
      whiteout = 1;
    }
    if (whiteout == 1) {
      a[i] = ' ';
    }
  }
}

/***=======================================================================***/
/*** ALPHANUMERIC: make sure a character is alphanumeric; if not, make it  ***/
/***               whitespace.                                             ***/
/***=======================================================================***/
void alphanumeric(char *a)
{
  if (isalphanum(*a) == 0) {
    *a = ' ';
  }
}

/***=======================================================================***/
/*** ISALPHANUM: return 0 if a character is not 'alphanumeric' by GRAPPLE  ***/
/***             standards.                                                ***/
/***=======================================================================***/
int isalphanum(char a)
{
  if ((a < 'a' || a > 'z') && (a < 'A' || a > 'Z') && (a < '0' || a > '9') &&
      (a != '*') && (a != '\'')) {
    return 0;
  }
  else {
    return 1;
  }
}

/***=======================================================================***/
/*** REALIGN_NAME: function for realigning an atom name of four characters ***/
/***               to look like the original PDB format (one space goes in ***/
/***               front, unless the first character is a letter or the    ***/
/***               name has four characters).                              ***/
/***=======================================================================***/
void realign_name(char* atom_name)
{
  int i;

  if (atom_name[0] != ' ' && atom_name[3] == ' ' &&
      (atom_name[0] < '0' || atom_name[0] > '9')) {
    for (i = 3; i > 0; i--) {
      atom_name[i] = atom_name[i-1];
    }
    atom_name[0] = ' ';
  }
}

/***=======================================================================***/
/*** CopyAtoms: copies atoms from one PDB to another.                      ***/
/***=======================================================================***/
void CopyAtoms(pdb *add, int astart, int aend, pdb *base, int bstart, int chn)
{
  int i, j, ib;

  for (i = astart; i < aend; i++) {
    ib = (i-astart) + bstart;
    for (j = 0; j < 4; j++) {
      base->atom_names[4*ib+j] = add->atom_names[4*i+j];
      base->res_names[4*ib+j] = add->res_names[4*i+j];
    }
    for (j = 0; j < 3; j++) {
      base->crds[3*ib+j] = add->crds[3*i+j];
    }
    base->atom_nums[ib] = add->atom_nums[i];
    base->res_nums[ib] = add->res_nums[i];
    base->chain[ib] = (chn > 0) ? chn : add->chain[i];
  }
}

/***=======================================================================***/
/*** DeleteResidue: delete a residue from a pdb structure, moving termini  ***/
/***                and residue limits as appropriate.                     ***/
/***=======================================================================***/
void DeleteResidue(pdb *mol, int r)
{
  int i, j, i3, i4, rlen, r3len, r4len;

  /*** Advance all atoms ***/
  rlen = mol->res_lims[r+1] - mol->res_lims[r];
  r4len = 4*rlen;
  r3len = 3*rlen;
  for (i = mol->res_lims[r+1]; i < mol->n_atoms; i++) {
    i3 = 3*i;
    i4 = 4*i;
    for (j = 0; j < 4; j++) {
      mol->atom_names[i4-r4len+j] = mol->atom_names[i4+j];
      mol->res_names[i4-r4len+j] = mol->res_names[i4+j];
    }
    for (j = 0; j < 3; j++) {
      mol->crds[i3-r3len+j] = mol->crds[i3+j];
    }
    mol->atom_nums[i-rlen] = mol->atom_nums[i];
    mol->res_nums[i-rlen] = mol->res_nums[i];
    mol->chain[i-rlen] = mol->chain[i];
  }

  /*** Update termini ***/
  for (i = 0; i < mol->n_ter; i++) {
    if (mol->termini[i] > mol->res_lims[r] &&
	mol->termini[i] <= mol->res_lims[r+1]) {
      for (j = i+1; j < mol->n_ter; j++) {
	mol->termini[j-1] = mol->termini[j];
      }
      mol->n_ter -= 1;
    }
    if (mol->termini[i] >= mol->res_lims[r+1]) {
      mol->termini[i] -= rlen;
    }
  }

  /*** Update residue limits ***/
  for (i = r+1; i < mol->n_res-1; i++) { 
    mol->res_lims[i] = mol->res_lims[i+1] - rlen;
  }

  /*** Update the numbers of atoms and residues ***/
  mol->n_atoms -= rlen;
  mol->n_res -= 1;
}

/***=======================================================================***/
/*** ModPdbRA: modify the numbers of residues and atoms so that they do    ***/
/***           not overload the fields.                                    ***/
/***=======================================================================***/
void ModPdbRA(pdb *w)
{
  int i, j;

  for (i = 0; i < w->n_res; i++) {
    for (j = w->res_lims[i]; j < w->res_lims[i+1]; j++) {
      w->res_nums[j] = i;
    }
  }
  for (i = 0; i < w->n_atoms; i++) {
    w->atom_nums[i] = i;
  }
  for (i = 0; i < w->n_atoms; i++) {
    w->atom_nums[i] = w->atom_nums[i] % 1000000;
    w->res_nums[i] = w->res_nums[i] % 10000;
  }
}

/***=======================================================================***/
/*** AugmentPDB: augment a PDB structure with additional information.      ***/
/***=======================================================================***/
void AugmentPDB(pdb *thispdb)
{
  int i, j, k;

  thispdb->aug = 1;
  thispdb->z = (int*)calloc(thispdb->n_atoms, sizeof(int));
  thispdb->mass = (double*)calloc(thispdb->n_atoms, sizeof(double));

  /*** Atomic numbers ***/
  for (i = 0; i < thispdb->n_atoms; i++) {
    j = 4*i;
    k = thispdb->atom_names[j];
    while(k > 47 && k < 58) {
      j++;
      k = thispdb->atom_names[j];
    }
    if (k == 72) {
      thispdb->z[i] = 1;
      thispdb->mass[i] = 1.0078;
    }
    if (k == 67) {
      thispdb->z[i] = 6;
      thispdb->mass[i] = 12.0000;
    }
    if (k == 79) {
      thispdb->z[i] = 8;
      thispdb->mass[i] = 15.9996;
    }
    if (k == 78) {
      thispdb->z[i] = 7;
      thispdb->mass[i] = 14.0067;
    }
    if (k == 83) {
      thispdb->z[i] = 16;
      thispdb->mass[i] = 32.0650;
    }
  }
}
