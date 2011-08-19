/*

COPYRIGHT 1995 Rutgers - The State University of New Jersey

This software is provided WITHOUT WARRANTY OF MERCHANTABILITY OR
FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR
IMPLIED.  RUTGERS MAKE NO REPRESENTATION OR WARRANTY THAT THE
SOFTWARE WILL NOT INFRINGE ANY PATENT, COPYRIGHT OR OTHER
PROPRIETARY RIGHT.

The user of this software shall indemnify, hold harmless and defend
Rutgers, its governors, trustees, officers, employees, students,
agents and the authors against any and all claims, suits,
losses, liabilities, damages, costs, fees, and expenses including
reasonable attorneys' fees resulting from or arising out of the
use of this software.  This indemnification shall include, but is
not limited to, any and all claims alleging products liability.

This software may be used only for not-for-profit educational
and research purposes.

*/

/* ****************************************************
 * range.h: Header for error handlers
******************************************************/
#include <stdlib.h>
#include <iostream.h>

#ifndef H_RANGE
#define H_RANGE

#ifdef NO_RANGE_CHECK
#define Check(i) i
#else
#define Check(i) CheckIndex(i)
#endif

#ifndef NO_RANGE_CHECK
extern unsigned DefaultRangeErrorHandler(char *type, unsigned i, unsigned sz);
extern unsigned (*HandleRangeError)(char *type, unsigned i, unsigned sz);
#endif

#endif

