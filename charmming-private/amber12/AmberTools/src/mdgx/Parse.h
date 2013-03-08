#ifndef ParseHeadings
#define ParseHeadings

#include "MatrixDS.h"
#include "TopologyDS.h"
#include "CrdManipDS.h"

int CountWords(char* line);

cmat ParseWords(char* line);

void RemoveWhiteSpace(char* a, int asize);

void EqualSpace(char* line);

void RemoveComments(char* line);

void NixCommaCarriage(char* line);

int AdvanceToSegment(FILE *inp, char* segname, int scan0);

int DetectNamelistEnd(char* line, char* errmsg);

void SeekString(cmat L, char* val, char* sname, char* salias);

cmat SeekNString(cmat L, cmat* val, int* fspec, char* sname, char* salias);

void SeekReal(cmat L, double *val, char* sname, char* salias);

void SeekNReal(cmat L, double* val, char* sname, char* salias, int maxidx);

void SeekInt(cmat L, int *val, char* sname, char* salias);

void SeekLLInt(cmat L, long long int *val, char* sname, char* salias);

int* ParseAmbMask(char* maskstr, prmtop *tp, coord *crd);

FILE* FOpenSafe(char* fname, int ovrwrt);

long long int ReadNumericalShorthand(char* numstr);

#endif
