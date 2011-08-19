/* ************************************************************** *
 *  array_construct.h: Header file containing inline routines
 *              to help support constructing elements of arrays.
 *
 *       Adapted and modified from placemnt.h 
 *          Practical Data Structures in C++
 *          Bryan Flamig, Azarona Software/John Wiley & Sons, Inc
 * ************************************************************** */
#include <sys/types.h>

#ifndef H_ARRAY_CONSTRUCT
#define H_ARRAY_CONSTRUCT

/* ****************************************************************
 *  Proto types of the inline functions 
 * ****************************************************************/

template<class TYPE> void CopyElement(TYPE *p, const TYPE &x);

#ifndef H_ARRAY_CONSTRUCT_P
#ifndef HAVE_PLACEMENT_NEW
inline void *operator new(size_t, void *p) 
/* ---------------------------------------------------------------
 *  This "placement new operator" will cause an object
 *  to be constructed at a specified address p.
 * ---------------------------------------------------------------*/
{ 
  return p; 
}
#endif /* HAVE_PLACEMENT_NEW not defined */
#endif

#ifdef INCL_TEMPLATE_SRC
#include "CifArrayConstruct.C"
#endif

#endif
