/*

COPYRIGHT 1995-1997 Rutgers - The State University of New Jersey

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

/* **************************************************************
   rangeerr.C: Error handler for template class array or string
 ****************************************************************/
#include "range.h"

#ifndef NO_RANGE_CHECK
unsigned (*HandleRangeError)(char *type, unsigned i, unsigned sz) = DefaultRangeErrorHandler;

unsigned DefaultRangeErrorHandler(char *type, unsigned i, unsigned sz)
{
  cout << " Subscript of " << type << ' ' << i << " out of range (0, " << (sz-1) << ")\n";
  exit(1);
  return 0; // Not used in this function
}
#endif
