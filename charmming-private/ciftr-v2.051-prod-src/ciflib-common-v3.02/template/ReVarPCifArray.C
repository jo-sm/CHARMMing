/* ***************************************************************
   ReVarPCifArray.c: Methods for resizable, variable-length arrays

      Adapted and modified from rvarray.mth 
         Practical Data Structures in C++
         Bryan Flamig, Azarona Software/John Wiley & Sons, Inc
 * ************************************************************** */

template<class TYPE> int ReVarPCifArray<TYPE>::Realloc(unsigned new_dimlen, int keep)
// ---------------------------------------------------------------
// Redimensions the array to the new dimension length new_dimlen.
// If grow_by = 0, or a can't get more memory, a 0 is returned, 
// else, a 1 is returned. If keep == 1, data is left intact 
// via copying of elements to new space, otherwise, the
// elements are destroyed by calling their destructors. 
// If the new length is shorter than the old current length
// of the array, the orphaned data is destroyed.
// NOTE: If error occurs, the data is left intact.
// ---------------------------------------------------------------
{
  if (grow_by == 0) return 0;
  char *new_data = (char *) new char[new_dimlen * sizeof(TYPE)];

  if (new_data == 0) return 0;
  this->dimlen = new_dimlen;
  if (keep) {
     if (new_dimlen < this->len) {
        this->len = new_dimlen;
     }
     // Copy old data into new space
     memmove(new_data, this->data, this->len*sizeof(TYPE));
  }
  else {
    this->len = 0;
  }
  delete[] (char *)(this->data); // Delete old storage place
  this->data = (TYPE *) new_data;
  return 1;
}

template<class TYPE> void ReVarPCifArray<TYPE>::CopyN(const TYPE *s, unsigned n)
// ---------------------------------------------------------------
// Copies the data from s, resizing this array if n is 
// larger than the dimensioned length. Note that the resize 
// may fail, in which case the array keeps it's old size,
// and not all of the data will be copied.
// ---------------------------------------------------------------
{
  if (n > this->dimlen) {

    if (n < this->dimlen + grow_by * 2)
      Grow(grow_by *2);
    else
      Grow(n-this->dimlen);
  }
  VarPCifArray<TYPE>::CopyN(s, n);
}

template<class TYPE>
unsigned ReVarPCifArray<TYPE>::InsertNAt(unsigned p, const TYPE *s, unsigned n)
// ---------------------------------------------------------------
// Inserts the data pointed to by s into the array. Up to n
// elements are inserted, (truncating if necesary), starting at
// position p, (counted by zero). The array will grow in size
// if needed and grow_by != 0. The size it will grow
// will be the larger of grow_by and the room_needed.
// If p >= len, then the data is concatenated on the end.
// Returns number of elements inserted, or 0 if error.
// ---------------------------------------------------------------
{
  unsigned room, needed;
  
  room = this->dimlen - this->len;
  if (n > room && grow_by != 0) {
     needed = n - room;

     if (needed < (unsigned int)(grow_by * 2))
       needed = grow_by * 2;

     if (needed < (unsigned int)(grow_by)) needed = grow_by;
     // Grow the array by needed amount. Note that the
     // growth may fail, in which case the array stays
     // the same size, with data intact.
     // s1 copying is added by Shuhsin on Dec. 1, 94
     // It is just in case s belong to one of the data in 
     // the array, and because of reallocation, s might 
      // get screwed up
     Grow(needed);
   return VarPCifArray<TYPE>::InsertNAt(p, s, n);
   }
  else
    return VarPCifArray<TYPE>::InsertNAt(p, s, n);
}


template<class TYPE> void ReVarPCifArray<TYPE>::DeleteElement() {
  delete[] (char *)(this->data);
};
