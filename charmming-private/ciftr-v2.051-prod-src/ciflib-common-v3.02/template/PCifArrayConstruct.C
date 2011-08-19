/* ************************************************************** *
 *  parray_construct.c: Header file containing inline routines
 *              to help support constructing elements of arrays.
 *
 *       Adapted and modified from placemnt.h 
 *          Practical Data Structures in C++
 *          Bryan Flamig, Azarona Software/John Wiley & Sons, Inc
 * ************************************************************** */

template<class TYPE> void CopyElem(TYPE *p, const TYPE &x)
/* ---------------------------------------------------------------
 *  This function constructs an element at place p. It
 *  assumes TYPE has a copy constructor. If TYPE doesn't
 *  have one but has a default constructor, you can
 *  override this function, and replace the body with
 *  a call to the default constructor fby an assignment:
 *  { new(p) TYPE; *p = x; }
 *  Or if TYPE is a built-in type, you may replace the
 *  body with simply: { *p = x; }.
 * ---------------------------------------------------------------*/
{
  *p = x;
}


