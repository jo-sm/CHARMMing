#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Constants.h"
#include "Random.h"

#include <sys/time.h>

/***=======================================================================***/
/*** SetPRNG: this sets the pseudo-random number generator based on a      ***/
/***          seed.  If the seed is less than zero, then the seed is taken ***/
/***          from the wall time and returned.                             ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   seed:  pointer to a long integer                                    ***/
/***=======================================================================***/
long int SetPRNG(int *seed)
{
  struct timeval trng;

  if (*seed < 0) {
    gettimeofday(&trng, NULL);
    *seed = 1.0e6*(trng.tv_sec % 1000) + trng.tv_usec;
  }

  return *seed;
}

/***=======================================================================***/
/*** RAN2: function for returning a single random number from a uniform    ***/
/***       distribution in the range (0, 1), exclusive of the endpoints.   ***/
/***       This code is taken from [REF]:                                  ***/
/***                                                                       ***/
/***       William H. Press, Saul A. Teukolsky, William T. Vetterling, and ***/
/***       Brian P. Flannery.  Numerical Recipes in C, Second Edition.     ***/
/***       Cambridge University Press, 1992.                               ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   idum:  pointer to a long unsigned integer                           ***/
/***=======================================================================***/
double ran2(long *idum) 
{ 
  int j;
  long k;
  static long idum2=123456789; 
  static long iy=0;
  static long iv[NTAB]; 
  double temp;

  if (*idum <= 0) { 
    if (-(*idum) < 1) *idum=1; 
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) { 
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1; 
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1; 
  *idum=IA1*(*idum-k*IQ1)-k*IR1; 
  if (*idum < 0) *idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; 
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2; 
  iv[j] = *idum;
  if (iy < 1) iy += IMM1; 
  if ((temp=AM*iy) > RNMX) return RNMX; 
  else return temp;
}

/***=======================================================================***/
/*** GaussBoxMuller: obtain a Gaussian distribution of random numbers from ***/
/***                 a uniform one.  The mean of the Gaussian is 0.0 and   ***/
/***                 its width is 1.0.  This function returns only one     ***/
/***                 number; although it can generate two they would be    ***/
/***                 correlated with each other.                           ***/
/***                                                                       ***/
/*** Arguments:                                                            ***/
/***   idum:  pointer to a long unsigned (typically 32-bit) integer        ***/
/***=======================================================================***/
double GaussBoxMuller(long *counter)
{
  double x1, x2, y;

  x1 = sqrt(-2.0*log(ran2(counter)));
  x2 = sin(TWOPI*ran2(counter));
  y = x1 * x2;

  return y;
}
