/* ***************************************************************
   VarCifArray.c: Variable-length array base class methods
  
        Adapted and modified from VarCifArray.mth 
           Practical Data Structures in C++
           Bryan Flamig, Azarona Software/John Wiley & Sons, Inc
  ***************************************************************/
#include "CifArrayConstruct.h"

template<class TYPE> void VarCifArray<TYPE>::CleanUpElements(unsigned s, unsigned f)
// ---------------------------------------------------------------
// Cleans up elements from start s to finish f-1, by calling
// their destructors explicitly. 
// ---------------------------------------------------------------
{
  TYPE *p = this->data+s;
  for (unsigned i=s; i<f; i++, p++)    p->~TYPE();
}

template<class TYPE> void VarCifArray<TYPE>::CopyN(const TYPE *s, unsigned slen)
// ---------------------------------------------------------------
// Copies as much data as possible from s into this
// array, truncating if need be. Copying is done via
// element to element assignment.
// ---------------------------------------------------------------
{
  unsigned i, tlen;

  // First, copy into the constructed portion.
  tlen = this->len;
  if (tlen > slen) {
     Shrink(tlen-slen); // We might be getting smaller
     tlen = slen;
  }
  for (i = 0; i<tlen; i++) (this->data)[i] = s[i];
  this->len = tlen;
  // Now, copy into the unconstructed portion via concatenation
  if (tlen < slen) Concatenate(s+i, slen-tlen);
}

template<class TYPE> unsigned VarCifArray<TYPE>::InsertNAt(unsigned p, const TYPE *s, unsigned n)
// ---------------------------------------------------------------
// Inserts the data pointed to by s into the array. Up to n
// elements are inserted, (truncating if necesary), starting at
// position p, (counted by zero). If p >= len, then the data
// is concatenated on the end.
// Return value: number of elements inserted.
// ---------------------------------------------------------------
{
  unsigned i;
  if (n > dimlen - this->len) { // Keep things in range
     n = dimlen - this->len;
  }
  if (p >= this->len) { // We're concatenating
     p = this->len;
  }
  else {
    // Make room for inserted data somewhere in the middle.
    memmove(this->data+p+n, this->data+p, (this->len-p)*sizeof(TYPE));
  }
  // Copy in source
  this->len += n;

  for(i = 0; i<n; i++) {
     // Copy-construct the element in place (it's presumed 
     // to be unconstructed at this time.)
    CopyElement(this->data+p+i, s[i]);
  }

  return n;
}

template<class TYPE> unsigned VarCifArray<TYPE>::DeleteAt(unsigned p, unsigned n)
// ---------------------------------------------------------------
// Deletes up to n elements from array at pos'n p.  If p
// is out of range nothing happens. If n == 0, nothing happens.
// Return value:  number of elements deleted.
// ---------------------------------------------------------------
{
  long pel; // long is used to prevent overflows during adds
  unsigned pe;

  if (p < this->len && n != 0) {
     // We may be chopping off the end, so keep in range
     pel = long(p) + long(n);
     if ((unsigned long) pel > this->len) pel = long(this->len);
     pe = unsigned(pel); // Guaranteed to be in range
     n = pe-p;
     // Now, clean up each element deleted
     CleanUpElements(p, pe);
     // Now, move elements up to take their place
     memmove(this->data+p, this->data+pe, (this->len-pe)*sizeof(TYPE));
     this->len -= n;
  } else n = 0;
  return n;
}

template<class TYPE> unsigned VarCifArray<TYPE>::Shrink(unsigned n)
// ---------------------------------------------------------------
// Shrink the array by n elements.
// Return value: number of elements had been shrunk.
// ---------------------------------------------------------------
{
  if (n > this->len) n = this->len;
  return DeleteAt(this->len-n, n);
}

#ifndef NO_RANGE_CHECK

template<class TYPE> unsigned VarCifArray<TYPE>::CheckIndex(unsigned int i) const
// ---------------------------------------------------------------
// Check for index being in bounds. If not in
// bounds, call the error handler.
// ---------------------------------------------------------------
{
  if (i >= dimlen) 
    i = HandleRangeErr("VarCifArray", i, dimlen);
  return i;
}

#endif

