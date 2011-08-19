/* ***************************************************************
   DymVarCifArray.h: Variable-length dynamically allocated
                  array class template.
 
       Adapted and modified from dvarray.h 
          Practical Data Structures in C++
          Bryan Flamig, Azarona Software/John Wiley & Sons, Inc
 * ***************************************************************/
#ifndef H_DYMVARARRAY
#define H_DYMVARARRAY

#include "VarCifArray.h"

#define INLINE

template<class TYPE> class DymVarCifArray : public VarCifArray<TYPE> {
public:
  DymVarCifArray();
  DymVarCifArray(unsigned n);
  DymVarCifArray(const DymVarCifArray<TYPE> &s);
  DymVarCifArray(const CifArray<TYPE> &s);
  DymVarCifArray(const TYPE *s, unsigned n);
  virtual ~DymVarCifArray();
  DymVarCifArray<TYPE> &operator=(const CifArray<TYPE> &s);
  DymVarCifArray<TYPE> &operator=(const DymVarCifArray<TYPE> &s);
};

template<class TYPE> INLINE DymVarCifArray<TYPE>::DymVarCifArray()
: VarCifArray<TYPE>((TYPE *)new char[8*sizeof(TYPE)], 8) 
// ---------------------------------------------------------------
// Allocates for a variable length array of n max elements.
// ---------------------------------------------------------------
{
  memset(this, 0, 8*sizeof(TYPE));
}
template<class TYPE> INLINE DymVarCifArray<TYPE>::DymVarCifArray(unsigned n)
: VarCifArray<TYPE>((TYPE *)new char[n*sizeof(TYPE)], n) 
// ---------------------------------------------------------------
// Allocates for a variable length array of n max elements.
// ---------------------------------------------------------------
{ 
}

template<class TYPE>
INLINE DymVarCifArray<TYPE> &DymVarCifArray<TYPE>::operator=(const CifArray<TYPE> &s)
// ---------------------------------------------------------------
// Note that source of assignment can be any array type.
// Traps assignment to self.
// ---------------------------------------------------------------
{ 
  if (this != &s) Copy(s); 
  return *this; 
}

template<class TYPE>
INLINE DymVarCifArray<TYPE> &DymVarCifArray<TYPE>::operator=(const DymVarCifArray<TYPE> &s)
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
#include "DymVarCifArray.C"
#endif

#endif


