/* ************************************************************* *
   PCifArray.c: Basic array class template implementation.
  
        Adapted and modified from array.mth 
           Practical Data Structures in C++
           Bryan Flamig, Azarona Software/John Wiley & Sons, Inc
 * *************************************************************** */
#include <stdio.h>
#include "range.h"

#ifndef NO_RANGE_CHECK

template<class TYPE> unsigned PCifArray<TYPE>::CheckIndex(unsigned int i) const
// ---------------------------------------------------------------
//  Check for index being in bounds. If not in
//  bounds, call the error handler.
//  Return value: index;
// ---------------------------------------------------------------
{
  if (i >= DimLength()) 
    i = HandleRangeError("PCifArray", i, DimLength());

  return i;
}

#endif

template<class TYPE>  void PCifArray<TYPE>::CopyN(const TYPE *s, unsigned slen)
// ---------------------------------------------------------------
//  Copies as much data as possible from s into this 
//  array, truncating if need be. Element by element
//  assignment is used for the copying.
// ---------------------------------------------------------------
{
  unsigned clen = DimLength();
  if (clen > slen) clen = slen;
  for (unsigned i = 0; i<clen; i++) data[i] = s[i];
}








