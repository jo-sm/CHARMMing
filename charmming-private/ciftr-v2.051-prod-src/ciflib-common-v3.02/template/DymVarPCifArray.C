/* ***************************************************************
   DymVarPCifArray.c: Variable-length dynamically allocated
                  array class methods.
 
       Adapted and modified from rvarray.mth 
          Practical Data Structures in C++
          Bryan Flamig, Azarona Software/John Wiley & Sons, Inc
 * ***************************************************************/



template<class TYPE> DymVarPCifArray<TYPE>::DymVarPCifArray(const DymVarPCifArray<TYPE> &s)
: VarPCifArray<TYPE>((TYPE *)new char[s.DimLength()*sizeof(TYPE)], 
// ---------------------------------------------------------------
// Copy constructor. Note that we first allocate an
// array of unconstructed bits. Using concatenate
// here effectively does the element-by-element copy
// construction for us.
// ---------------------------------------------------------------
               s.DimLength())
{
  memset(this->data, 0, s.DimLength()* sizeof(TYPE));
  Concatenate(s);
}

template<class TYPE>
DymVarPCifArray<TYPE>::DymVarPCifArray(const PCifArray<TYPE> &s)
: VarPCifArray<TYPE>((TYPE *)new char[s.DimLength()*sizeof(TYPE)], 
               s.DimLength())
// ---------------------------------------------------------------
// Constructor to make a dynamically allocated, variable
// length copy of another type of array. See comments to
// copy constructor.
// ---------------------------------------------------------------
{
  memset(this->data, 0, s.DimLength()* sizeof(TYPE));
  Concatenate(s);
}

template<class TYPE>
DymVarPCifArray<TYPE>::DymVarPCifArray(const TYPE *s, unsigned n)
// ---------------------------------------------------------------
// Constructor to make a dynamically allocated, variable
// length copy of some low-level memory. See comments to
// copy constructor.
// ---------------------------------------------------------------
: VarPCifArray<TYPE>((TYPE *)new char[n*sizeof(TYPE)], n)
{
  memset(this->data, 0, n*sizeof(TYPE));
  Concatenate(s, n);
}

template<class TYPE> DymVarPCifArray<TYPE>::~DymVarPCifArray()
// ---------------------------------------------------------------
//  Destructor
// ---------------------------------------------------------------
{
  // Delete the storage
  delete[] (char *)(this->data);
}
