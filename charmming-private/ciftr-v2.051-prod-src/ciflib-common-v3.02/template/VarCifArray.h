/* ***************************************************************
   VarCifArray.h: Variable-length array base class template
  
        Adapted and modified from VarCifArray.h
           Practical Data Structures in C++
           Bryan Flamig, Azarona Software/John Wiley & Sons, Inc
 * ***************************************************************/
#ifndef H_VARARRAY
#define H_VARARRAY
#include <string.h>
#include "CifArray.h"

#define INLINE

// ***************************************************************
//  WARNING: In the following templates, TYPE is assumed
//  to have a copy constructor. If it doesn't, see the
//  function CopyElement() in array_construct.h for what to do.
//
//  ALSO: If the types you are going to use have no
//  destructor/DeleteElement() function, you may optionally 
//  override the template function DeleteElement() in 
//  array_construct.h accordingly. With some compilers, 
//  you may be required to do this.
//
//  NOTE: Element-to-element assignments are used
//  for all inter-array copying.
//
// 
//  Variable-length array base class template. These arrays
//  have data allocated elsewhere.
// 
// ***************************************************************/

template<class TYPE> class VarCifArray : public CifArray<TYPE> {
 protected:
  unsigned dimlen;

 public:
  VarCifArray(TYPE *m, unsigned dim_sz);
  void CleanUpElements(unsigned s, unsigned f);
  virtual void CopyN(const TYPE *s, unsigned slen);
  VarCifArray<TYPE> &operator=(const CifArray<TYPE> &s);
  VarCifArray<TYPE> &operator=(const VarCifArray<TYPE> &s);
#ifndef NO_RANGE_CHECK
  unsigned CheckIndex(unsigned i) const;
#endif
  virtual  unsigned InsertNAt(unsigned p, const TYPE *s, unsigned n);
  unsigned Shrink(unsigned n);

  virtual unsigned DimLength() const { return dimlen; }
  unsigned InsertAt(unsigned p, const TYPE &s) { return InsertNAt(p, &s, 1); }
  unsigned InsertAt(unsigned p, const CifArray<TYPE> &a) 
    { return InsertNAt(p, a.Data(), a.Length()); }
  unsigned Concatenate(const TYPE *s, unsigned n=1) { return InsertNAt(this->len, s, n); }
  unsigned Concatenate(const CifArray<TYPE> &a)   { return InsertAt(this->len, a); }
  unsigned Add(const TYPE &s) { return InsertNAt(this->len, &s, 1);}
  void Clear() { DeleteAt(0, this->len); }

  unsigned DeleteAt(unsigned p, unsigned n=1);
};


template<class TYPE> 
INLINE VarCifArray<TYPE>::VarCifArray(TYPE *m, unsigned dim_sz): CifArray<TYPE>(m, 0) 
// ---------------------------------------------------------------
// Constructor that wraps a variable length array around
// low-level memory m, which has room for dim_sz elements.
// The memory is presumed to be unconstructed bits.
// ---------------------------------------------------------------
{ 
  dimlen = dim_sz; 
}

template<class TYPE> 
INLINE VarCifArray<TYPE> &VarCifArray<TYPE>::operator=(const CifArray<TYPE> &s)
// ---------------------------------------------------------------
// Note that source of assignment can be any array type.
// Traps assignment to self.
// ---------------------------------------------------------------
{ 
  if (this != &s) Copy(s); 
  return *this; 
}

template<class TYPE> 
INLINE VarCifArray<TYPE> &VarCifArray<TYPE>::operator=(const VarCifArray<TYPE> &s)
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
#include "VarCifArray.C"
#endif


#endif



