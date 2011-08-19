/* ***************************************************************
   ReVarPCifArray.h: Resizable, variable-length array class template
 
       Adapted and modified from rvarray.h 
          Practical Data Structures in C++
          Bryan Flamig, Azarona Software/John Wiley & Sons, Inc
 * ***************************************************************/

#ifndef H_REVARARRAY_P
#define H_REVARARRAY_P
#include "DymVarPCifArray.h"

#define INLINE

template<class TYPE> class ReVarPCifArray : public DymVarPCifArray<TYPE>{
 protected:
  int grow_by;

 public:
  ReVarPCifArray();
  ReVarPCifArray(unsigned n, int gb);
  ReVarPCifArray(const ReVarPCifArray<TYPE> &s);
  ReVarPCifArray(const PCifArray<TYPE> &s);
  ReVarPCifArray(const TYPE *s, unsigned n);
  virtual void CopyN(const TYPE *s, unsigned n);
  ReVarPCifArray<TYPE> &operator=(const PCifArray<TYPE> &s);
  ReVarPCifArray<TYPE> &operator=(const ReVarPCifArray<TYPE> &s);
  virtual unsigned InsertNAt(unsigned p, const TYPE *s, unsigned n=1);
  virtual unsigned Add(const TYPE s);
  int Realloc(unsigned new_dimlen, int keep=1);
  int GrowBy(unsigned amt);
  int Grow();
  int Grow(unsigned grow_by_amount);
  int FreeExtra();
  void ChgGrowByInc(int gb);
  void DeleteElement();
};

template<class TYPE> INLINE ReVarPCifArray<TYPE>::ReVarPCifArray()
// ---------------------------------------------------------------
// Constructs a resizable array with a grow_by increment of gb.
// ---------------------------------------------------------------
:  DymVarPCifArray<TYPE>(1)
{ 
  grow_by = 1;
}
template<class TYPE> INLINE ReVarPCifArray<TYPE>::ReVarPCifArray(unsigned n, int gb)
// ---------------------------------------------------------------
// Constructs a resizable array with a grow_by increment of gb.
// ---------------------------------------------------------------
: DymVarPCifArray<TYPE>(n)
{ 
  grow_by = gb; 
}

template<class TYPE> INLINE ReVarPCifArray<TYPE>::ReVarPCifArray(const ReVarPCifArray<TYPE> &s)
// ---------------------------------------------------------------
// Copy constructor.
// ---------------------------------------------------------------
: DymVarPCifArray<TYPE>(s)
{  
  grow_by = s.grow_by; // Or should we say 0?
}

template<class TYPE> INLINE ReVarPCifArray<TYPE>::ReVarPCifArray(const PCifArray<TYPE> &s)
// ---------------------------------------------------------------
// Constructor that makes a copy of the data in other
// types of array. Defaults to not allowing any resizing.
// ---------------------------------------------------------------
: DymVarPCifArray<TYPE>(s)
{  
  grow_by = 0; 
}

template<class TYPE>
INLINE ReVarPCifArray<TYPE>::ReVarPCifArray(const TYPE *s, unsigned n)
// ---------------------------------------------------------------
// Constructor that creates resizable array and copies the
// data from the low-level C array. Defaults to not allowing
// any resizing.
// ---------------------------------------------------------------
: DymVarPCifArray<TYPE>(s, n)
{ 
  grow_by = 0; 
}

template<class TYPE>
INLINE ReVarPCifArray<TYPE> &ReVarPCifArray<TYPE>::operator=(const PCifArray<TYPE> &s)
// ---------------------------------------------------------------
// Note that source of assignment can be any array type.
// Traps assignment to self.
// ---------------------------------------------------------------
{ 
  if (this != &s) Copy(s); 
  return *this; 
}

template<class TYPE>
INLINE ReVarPCifArray<TYPE> &ReVarPCifArray<TYPE>::operator=(const ReVarPCifArray<TYPE> &s)
// ---------------------------------------------------------------
// Note that source of assignment can be any array type.
// Traps assignment to self.
// ---------------------------------------------------------------
{ 
  if (this != &s) Copy(s); 
  return *this; 
}
template<class TYPE>
INLINE unsigned ReVarPCifArray<TYPE>::Add(const TYPE s) 
// ---------------------------------------------------------------
// Add an element onto the end of the arary.
// ---------------------------------------------------------------
{ 
  return InsertNAt(this->len, &s, 1); 
}

template<class TYPE>
INLINE int ReVarPCifArray<TYPE>::GrowBy(unsigned amt)
// ---------------------------------------------------------------
// Grows the array by the specifed amt.
// ---------------------------------------------------------------
{
  return Realloc(this->dimlen+amt);
}

template<class TYPE>
INLINE int ReVarPCifArray<TYPE>::Grow()
// ---------------------------------------------------------------
// Grows the array by the default amount.
// ---------------------------------------------------------------
{
  grow_by *= 2;
  return Realloc(this->dimlen+grow_by); 
}

template<class TYPE>
INLINE int ReVarPCifArray<TYPE>::Grow(unsigned grow_by_amount)
// ---------------------------------------------------------------
// Grows the array by the default amount.
// ---------------------------------------------------------------
{
  ChgGrowByInc(grow_by_amount);
  return Realloc(this->dimlen+grow_by); 
}

template<class TYPE>
INLINE int ReVarPCifArray<TYPE>::FreeExtra()
// ---------------------------------------------------------------
// Frees all unused space in the array.
// ---------------------------------------------------------------
{ 
  return Realloc(this->len);
}

template<class TYPE> INLINE void ReVarPCifArray<TYPE>::ChgGrowByInc(int gb)
{ 
  grow_by = gb; 
}
#undef INLINE

#ifdef INCL_TEMPLATE_SRC
#include "ReVarPCifArray.C"
#endif

#endif
