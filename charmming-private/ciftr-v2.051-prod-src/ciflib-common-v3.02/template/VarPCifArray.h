/* ***************************************************************
   VarPCifArray.h: Variable-length array base class template
  
        Adapted and modified from VarPCifArray.h
           Practical Data Structures in C++
           Bryan Flamig, Azarona Software/John Wiley & Sons, Inc
 * ***************************************************************/
#ifndef H_VARARRAY_P
#define H_VARARRAY_P
#include <string.h>
#include "PCifArray.h"

#define INLINE

// ***************************************************************
//  WARNING: In the following templates, TYPE is assumed
//  to have a copy constructor. If it doesn't, see the
//  function CopyElem() in array_construct.h for what to do.
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

template<class TYPE> class VarPCifArray : public PCifArray<TYPE> {
 protected:
  unsigned dimlen;

 public:
  VarPCifArray(TYPE *m, unsigned dim_sz);
  virtual void CopyN(const TYPE *s, unsigned slen);
  VarPCifArray<TYPE> &operator=(const PCifArray<TYPE> &s);
  VarPCifArray<TYPE> &operator=(const VarPCifArray<TYPE> &s);
#ifndef NO_RANGE_CHECK
  unsigned CheckIndex(unsigned i) const;
#endif
  virtual  unsigned InsertNAt(unsigned p, const TYPE *s, unsigned n);
  unsigned Shrink(unsigned n);

  virtual unsigned DimLength() const { return dimlen; }
  unsigned InsertAt(unsigned p, const TYPE &s) { return InsertNAt(p, &s, 1); }
  unsigned InsertAt(unsigned p, const PCifArray<TYPE> &a) 
    { return InsertNAt(p, a.Data(), a.Length()); }
  unsigned Concatenate(const TYPE *s, unsigned n=1) { return InsertNAt(this->len, s, n); }
  unsigned Concatenate(const PCifArray<TYPE> &a)   { return InsertAt(this->len, a); }
  unsigned Add(const TYPE s) { return InsertNAt(this->len, &s, 1);}
  void Clear() { DeleteAt(0, this->len); }

  unsigned DeleteAt(unsigned p, unsigned n=1);
};


template<class TYPE> 
INLINE VarPCifArray<TYPE>::VarPCifArray(TYPE *m, unsigned dim_sz): PCifArray<TYPE>(m, 0) 
// ---------------------------------------------------------------
// Constructor that wraps a variable length array around
// low-level memory m, which has room for dim_sz elements.
// The memory is presumed to be unconstructed bits.
// ---------------------------------------------------------------
{ 
  dimlen = dim_sz; 
}

template<class TYPE> 
INLINE VarPCifArray<TYPE> &VarPCifArray<TYPE>::operator=(const PCifArray<TYPE> &s)
// ---------------------------------------------------------------
// Note that source of assignment can be any array type.
// Traps assignment to self.
// ---------------------------------------------------------------
{ 
  if (this != &s) Copy(s); 
  return *this; 
}

template<class TYPE> 
INLINE VarPCifArray<TYPE> &VarPCifArray<TYPE>::operator=(const VarPCifArray<TYPE> &s)
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
#include "VarPCifArray.C"
#endif


#endif



