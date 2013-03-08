#ifndef TimingsHeaders
#define TimingsHeaders

#include "TimingsDS.h"
#include "pmeRecipDS.h"

void InitExecon(execon *tm);

void mdgxStartTimer(execon *tm);

double mdgxStopTimer(execon *tm);

void PrintTimingData(execon *tm, reccon *rcinp, FILE *outp);

#endif
