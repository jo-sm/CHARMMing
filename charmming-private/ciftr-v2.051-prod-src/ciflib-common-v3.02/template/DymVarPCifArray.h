/* ***************************************************************
   DymVarPCifArray.h: Variable-length dynamically allocated
                  array class template.
 
       Adapted and modified from dvarray.h 
          Practical Data Structures in C++
          Bryan Flamig, Azarona Software/John Wiley & Sons, Inc
 * ***************************************************************/
#ifndef H_DYMVARARRAY_P
#define H_DYMVARARRAY_P

#include "VarPCifArray.h"

#define INLINE

template<class TYPE> class DymVarPCifArray : public VarPCifArray<TYPE> {
public:
  DymVarPCifArray();
  DymVarPCifArray(unsigned n);
  DymVarPCifArray(const DymVarPCifArray<TYPE> &s);
  DymVarPCifArray(const PCifArray<TYPE> &s);
  DymVarPCifArray(const TYPE *s, unsigned n);
  virtual ~DymVarPCifArray();
  DymVarPCifArray<TYPE> &operator=(const PCifArray<TYPE> &s);
  DymVarPCifArray<TYPE> &operator=(const DymVarPCifArray<TYPE> &s);
};

template<class TYPE> INLINE DymVarPCifArray<TYPE>::DymVarPCifArray()
: VarPCifArray<TYPE>((TYPE *)new char[8*sizeof(TYPE)], 8) 
// ---------------------------------------------------------------
// Allocates for a variable length array of n max elements.
// ---------------------------------------------------------------
{
  memset(this, 0, 8*sizeof(TYPE));
}

template<class TYPE> INLINE DymVarPCifArray<TYPE>::DymVarPCifArray(unsigned n)
: VarPCifArray<TYPE>((TYPE *)new char[n*sizeof(TYPE)], n) 
// ---------------------------------------------------------------
// Allocates for a variable length array of n max elements.
// ---------------------------------------------------------------
{ 
}

template<class TYPE>
INLINE DymVarPCifArray<TYPE> &DymVarPCifArray<TYPE>::operator=(const PCifArray<TYPE> &s)
// ---------------------------------------------------------------
// Note that source of assignment can be any array type.
// Traps assignment to self.
// ---------------------------------------------------------------
{ 
  if (this != &s) Copy(s); 
  return *this; 
}

template<class TYPE>
INLINE DymVarPCifArray<TYPE> &DymVarPCifArray<TYPE>::operator=(const DymVarPCifArray<TYPE> &s)
// ---------------------------------------------------------------
// Note that source of assignment can be any array type.
// Traps assignment to self.
// ---------------------------------------------------------------
{ 
  if (this != &s) Copy(s); 
  return *this; 
}

#undef INLINE


#ifdef INCL_TEMPLATE_SRC
#include "DymVarPCifArray.C"
#endif

#endif


