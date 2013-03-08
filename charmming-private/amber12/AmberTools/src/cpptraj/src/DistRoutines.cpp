#include <cmath> //floor
#include <cstddef> // NULL
#include "DistRoutines.h"
#include "Constants.h" // DEGRAD

/* Frame::ClosestImage()
 * Given two coordinates A and B, determine the unit XYZ vector that points 
 * towards the closest image of B to A.
 * It is assumed the coordinates are already relative to the box center.
 */
/*
void Frame::ClosestImage(double *A, double *B, int *ixyz) {
  double halfBox[3];//, delta;
  int vectorA[3], vectorB[3], i;

  mprintf("DEBUG: CoordA     = %lf %lf %lf\n",A[0],A[1],A[2]);
  mprintf("DEBUG: CoordB     = %lf %lf %lf\n",B[0],B[1],B[2]);

  halfBox[0] = box[0] / 6.0; 
  halfBox[1] = box[1] / 6.0; 
  halfBox[2] = box[2] / 6.0;
  mprintf("DEBUG: Half box = %lf %lf %lf\n",halfBox[0],halfBox[1],halfBox[2]);

  // Vector A
  vectorA[0] = 0; vectorA[1] = 0; vectorA[2] = 0;
  for (i=0; i<3; i++) {
    if (A[i] < -halfBox[i]) vectorA[i] = -1;
    if (A[i] >  halfBox[i]) vectorA[i] =  1;
//    delta = A[i] - boxCenter[i];
//  if (delta > 0.0) vectorA[i] = 1;
//  if (delta < 0.0) vectorA[i] = -1;
//
  }
  mprintf("DEBUG:  VectorA = %2i %2i %2i\n",vectorA[0],vectorA[1],vectorA[2]);

  // NOT Vector B
  vectorB[0] = 0; vectorB[1] = 0; vectorB[2] = 0;
  for (i=0; i<3; i++) {
    if (B[i] < -halfBox[i]) vectorB[i] =  1; // NOT
    if (B[i] >  halfBox[i]) vectorB[i] = -1; // NOT
//    delta = B[i] - boxCenter[i];
//  if (delta > 0.0) vectorB[i] = -1; // NOT
//  if (delta < 0.0) vectorB[i] = 1;  // NOT
//
  }
  mprintf("DEBUG: !VectorB = %2i %2i %2i\n",vectorB[0],vectorB[1],vectorB[2]);

  // A & !B
  ixyz[0]=0; ixyz[1]=0; ixyz[2]=0;
  for (i=0; i<3; i++) {
    if (vectorA[i] == vectorB[i]) ixyz[i] = vectorA[i];
    //ixyz[i] = vectorA[i] & vectorB[i];
  }
}
*/

// MinImageNonOrtho2()
/** Given two sets of coordinates and reciprocal space information based on
  * the current non-orthorhombic box, return the shortest imaged distance^2
  * between the coordinates.
  * The integer coefficients describing the closest reflection in reciprocal
  * space will be placed in ixyz.
  */
double MinImageNonOrtho2(double *Coord1, double *Coord2, double *box, int origin, int *ixyz,
                         double *ucell, double *recip) {
  double min, f[3], f2[3];

  min = 100.0 * (box[0]*box[0]+box[1]*box[1]+box[2]*box[2]);

  //if (prnlev > 6) {
  //  fprintf(stdout, "ATOM      0  XXX A1      1     %7.3f %7.3f %7.3f\n",
  //          x1, y1, z1);
  //  fprintf(stdout, "ATOM      1  XXX A2      1     %7.3f %7.3f %7.3f\n",
  //          x2, y2, z2);
  //}

  f[0] = Coord1[0]*recip[0] + Coord1[1]*recip[1] + Coord1[2]*recip[2];
  f[1] = Coord1[0]*recip[3] + Coord1[1]*recip[4] + Coord1[2]*recip[5];
  f[2] = Coord1[0]*recip[6] + Coord1[1]*recip[7] + Coord1[2]*recip[8];

  f2[0] = Coord2[0]*recip[0] + Coord2[1]*recip[1] + Coord2[2]*recip[2];
  f2[1] = Coord2[0]*recip[3] + Coord2[1]*recip[4] + Coord2[2]*recip[5];
  f2[2] = Coord2[0]*recip[6] + Coord2[1]*recip[7] + Coord2[2]*recip[8];

  if (origin) {
    f[0] += 0.5;
    f[1] += 0.5;
    f[2] += 0.5;
    f2[0] += 0.5;
    f2[1] += 0.5;
    f2[2] += 0.5;
  }

  min = DIST2_ImageNonOrthoRecip(f, f2, min, ixyz, ucell);

  return min;
}

// Frame::DIST2_ImageNonOrtho()
/** Given two coordinates and reciprocal space information based on 
  * the current non-orthorhombic box, return the shortest imaged distance^2 
  * between the coordinates.
  */
double DIST2_ImageNonOrtho(double *a1, double *a2, double *ucell, double *recip) { 
// double closest2
  double f[3], f2[3];
  int ixyz[3];

  f[0] = a2[0]*recip[0] + a2[1]*recip[1] + a2[2]*recip[2];
  f[1] = a2[0]*recip[3] + a2[1]*recip[4] + a2[2]*recip[5];
  f[2] = a2[0]*recip[6] + a2[1]*recip[7] + a2[2]*recip[8];

  f2[0] = a1[0]*recip[0] + a1[1]*recip[1] + a1[2]*recip[2];
  f2[1] = a1[0]*recip[3] + a1[1]*recip[4] + a1[2]*recip[5];
  f2[2] = a1[0]*recip[6] + a1[1]*recip[7] + a1[2]*recip[8];

  return DIST2_ImageNonOrthoRecip(f, f2, -1.0, ixyz, ucell);
}

// DIST2_ImageNonOrthoRecip()
/** Given two coordinate sets in reciprocal space, return the minimum imaged
  * distance^2 between them.
  * If minIn is > 0.0 it is considered a possible minimum distance.
  * The integer coefficients describing the closest reflection in reciprocal
  * space will be placed in ixyz.
  */
double DIST2_ImageNonOrthoRecip(double *f, double *f2, double minIn, int *ixyz, double *ucell) { 
  //double closest2
  double fx, fy, fz, f2x, f2y, f2z, X_factor, Y_factor, Z_factor;
  double fxm1, fxp1, fym1, fyp1, fzm1, fzp1;
  double x,y,z,D,min;

  /*  NON-ORTHORHOMBIC CASE: find shortest distance in periodic reference
   *  This is a brute force check requiring up to 26 distance evaluations.
   *  It has been adapted to be smarter by returning the first distance that
   *  is shorter than the minimum possible distance between images.
   */

  fx = f[0] - floor(f[0]);
  fy = f[1] - floor(f[1]);
  fz = f[2] - floor(f[2]); 
  
  f2x = f2[0] - floor(f2[0]);
  f2y = f2[1] - floor(f2[1]);
  f2z = f2[2] - floor(f2[2]);

  // Precompute some factors
  X_factor = (f2x*ucell[0] + f2y*ucell[3] + f2z*ucell[6]);
  Y_factor = (f2x*ucell[1] + f2y*ucell[4] + f2z*ucell[7]);
  Z_factor = (f2x*ucell[2] + f2y*ucell[5] + f2z*ucell[8]);

  fxm1 = fx - 1; fxp1 = fx + 1;
  fym1 = fy - 1; fyp1 = fy + 1;
  fzm1 = fz - 1; fzp1 = fz + 1;

  double fxm1u0 = fxm1 * ucell[0];
  double fxu0   = fx   * ucell[0];
  double fxp1u0 = fxp1 * ucell[0];
  double fxm1u1 = fxm1 * ucell[1];
  double fxu1   = fx   * ucell[1];
  double fxp1u1 = fxp1 * ucell[1];
  double fxm1u2 = fxm1 * ucell[2];
  double fxu2   = fx   * ucell[2];
  double fxp1u2 = fxp1 * ucell[2];

  double fym1u3 = fym1 * ucell[3];
  double fyu3   = fy   * ucell[3];
  double fyp1u3 = fyp1 * ucell[3];
  double fym1u4 = fym1 * ucell[4];
  double fyu4   = fy   * ucell[4];
  double fyp1u4 = fyp1 * ucell[4];
  double fym1u5 = fym1 * ucell[5];
  double fyu5   = fy   * ucell[5];
  double fyp1u5 = fyp1 * ucell[5];

  double fzm1u6 = fzm1 * ucell[6];
  double fzu6   = fz   * ucell[6];
  double fzp1u6 = fzp1 * ucell[6];
  double fzm1u7 = fzm1 * ucell[7];
  double fzu7   = fz   * ucell[7];
  double fzp1u7 = fzp1 * ucell[7];
  double fzm1u8 = fzm1 * ucell[8];
  double fzu8   = fz   * ucell[8];
  double fzp1u8 = fzp1 * ucell[8];

  // Calc ix iy iz = 0 case
  x = (fxu0 + fyu3 + fzu6) - X_factor;
  y = (fxu1 + fyu4 + fzu7) - Y_factor;
  z = (fxu2 + fyu5 + fzu8) - Z_factor;
  // DEBUG
  //fprintf(stdout,"DEBUG: a2: fx  fy  fz  = %lf %lf %lf\n",fx,fy,fz);
  //fprintf(stdout,"DEBUG: a1: f2x f2y f2z = %lf %lf %lf\n",f2x,f2y,f2z);
  min = (x*x) + (y*y) + (z*z);

  if (minIn > 0.0 && minIn < min) min = minIn;

  ixyz[0] = 0;
  ixyz[1] = 0;
  ixyz[2] = 0;

  // -1 -1 -1
  x = (fxm1u0 + fym1u3 + fzm1u6) - X_factor;
  y = (fxm1u1 + fym1u4 + fzm1u7) - Y_factor;
  z = (fxm1u2 + fym1u5 + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] = -1; ixyz[2] = -1; }
  // -1 -1  0
  x = (fxm1u0 + fym1u3 + fzu6  ) - X_factor;
  y = (fxm1u1 + fym1u4 + fzu7  ) - Y_factor;
  z = (fxm1u2 + fym1u5 + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] = -1; ixyz[2] =  0; }
  // -1 -1 +1
  x = (fxm1u0 + fym1u3 + fzp1u6) - X_factor;
  y = (fxm1u1 + fym1u4 + fzp1u7) - Y_factor;
  z = (fxm1u2 + fym1u5 + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] = -1; ixyz[2] =  1; }
  // -1  0 -1
  x = (fxm1u0 + fyu3   + fzm1u6) - X_factor;
  y = (fxm1u1 + fyu4   + fzm1u7) - Y_factor;
  z = (fxm1u2 + fyu5   + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] =  0; ixyz[2] = -1; }
  // -1  0  0
  x = (fxm1u0 + fyu3   + fzu6  ) - X_factor;
  y = (fxm1u1 + fyu4   + fzu7  ) - Y_factor;
  z = (fxm1u2 + fyu5   + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] =  0; ixyz[2] =  0; }
  // -1  0 +1
  x = (fxm1u0 + fyu3   + fzp1u6) - X_factor;
  y = (fxm1u1 + fyu4   + fzp1u7) - Y_factor;
  z = (fxm1u2 + fyu5   + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] =  0; ixyz[2] =  1; }
  // -1 +1 -1
  x = (fxm1u0 + fyp1u3 + fzm1u6) - X_factor;
  y = (fxm1u1 + fyp1u4 + fzm1u7) - Y_factor;
  z = (fxm1u2 + fyp1u5 + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] =  1; ixyz[2] = -1; }
  // -1 +1  0
  x = (fxm1u0 + fyp1u3 + fzu6  ) - X_factor;
  y = (fxm1u1 + fyp1u4 + fzu7  ) - Y_factor;
  z = (fxm1u2 + fyp1u5 + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] =  1; ixyz[2] =  0; }
  // -1 +1 +1
  x = (fxm1u0 + fyp1u3 + fzp1u6) - X_factor;
  y = (fxm1u1 + fyp1u4 + fzp1u7) - Y_factor;
  z = (fxm1u2 + fyp1u5 + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] =  1; ixyz[2] =  1; }

  //  0 -1 -1
  x = (fxu0   + fym1u3 + fzm1u6) - X_factor;
  y = (fxu1   + fym1u4 + fzm1u7) - Y_factor;
  z = (fxu2   + fym1u5 + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] = -1; ixyz[2] = -1; }
  //  0 -1  0
  x = (fxu0   + fym1u3 + fzu6  ) - X_factor;
  y = (fxu1   + fym1u4 + fzu7  ) - Y_factor;
  z = (fxu2   + fym1u5 + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] = -1; ixyz[2] =  0; }
  //  0 -1 +1
  x = (fxu0   + fym1u3 + fzp1u6) - X_factor;
  y = (fxu1   + fym1u4 + fzp1u7) - Y_factor;
  z = (fxu2   + fym1u5 + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] = -1; ixyz[2] =  1; }
  //  0  0 -1
  x = (fxu0   + fyu3   + fzm1u6) - X_factor;
  y = (fxu1   + fyu4   + fzm1u7) - Y_factor;
  z = (fxu2   + fyu5   + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] =  0; ixyz[2] = -1; }
  //  0  0  0
  //  0  0 +1
  x = (fxu0   + fyu3   + fzp1u6) - X_factor;
  y = (fxu1   + fyu4   + fzp1u7) - Y_factor;
  z = (fxu2   + fyu5   + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] =  0; ixyz[2] =  1; }
  //  0 +1 -1
  x = (fxu0   + fyp1u3 + fzm1u6) - X_factor;
  y = (fxu1   + fyp1u4 + fzm1u7) - Y_factor;
  z = (fxu2   + fyp1u5 + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] =  1; ixyz[2] = -1; }
  //  0 +1  0
  x = (fxu0   + fyp1u3 + fzu6  ) - X_factor;
  y = (fxu1   + fyp1u4 + fzu7  ) - Y_factor;
  z = (fxu2   + fyp1u5 + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] =  1; ixyz[2] =  0; }
  //  0 +1 +1
  x = (fxu0   + fyp1u3 + fzp1u6) - X_factor;
  y = (fxu1   + fyp1u4 + fzp1u7) - Y_factor;
  z = (fxu2   + fyp1u5 + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] =  1; ixyz[2] =  1; }

  // +1 -1 -1
  x = (fxp1u0 + fym1u3 + fzm1u6) - X_factor;
  y = (fxp1u1 + fym1u4 + fzm1u7) - Y_factor;
  z = (fxp1u2 + fym1u5 + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] = -1; ixyz[2] = -1; }
  // +1 -1  0
  x = (fxp1u0 + fym1u3 + fzu6  ) - X_factor;
  y = (fxp1u1 + fym1u4 + fzu7  ) - Y_factor;
  z = (fxp1u2 + fym1u5 + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] = -1; ixyz[2] =  0; }
  // +1 -1 +1
  x = (fxp1u0 + fym1u3 + fzp1u6) - X_factor;
  y = (fxp1u1 + fym1u4 + fzp1u7) - Y_factor;
  z = (fxp1u2 + fym1u5 + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] = -1; ixyz[2] =  1; }
  // +1  0 -1
  x = (fxp1u0 + fyu3   + fzm1u6) - X_factor;
  y = (fxp1u1 + fyu4   + fzm1u7) - Y_factor;
  z = (fxp1u2 + fyu5   + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] =  0; ixyz[2] = -1; }
  // +1  0  0
  x = (fxp1u0 + fyu3   + fzu6  ) - X_factor;
  y = (fxp1u1 + fyu4   + fzu7  ) - Y_factor;
  z = (fxp1u2 + fyu5   + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] =  0; ixyz[2] =  0; }
  // +1  0 +1
  x = (fxp1u0 + fyu3   + fzp1u6) - X_factor;
  y = (fxp1u1 + fyu4   + fzp1u7) - Y_factor;
  z = (fxp1u2 + fyu5   + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] =  0; ixyz[2] =  1; }
  // +1 +1 -1
  x = (fxp1u0 + fyp1u3 + fzm1u6) - X_factor;
  y = (fxp1u1 + fyp1u4 + fzm1u7) - Y_factor;
  z = (fxp1u2 + fyp1u5 + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] =  1; ixyz[2] = -1; }
  // +1 +1  0
  x = (fxp1u0 + fyp1u3 + fzu6  ) - X_factor;
  y = (fxp1u1 + fyp1u4 + fzu7  ) - Y_factor;
  z = (fxp1u2 + fyp1u5 + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] =  1; ixyz[2] =  0; }
  // +1 +1 +1
  x = (fxp1u0 + fyp1u3 + fzp1u6) - X_factor;
  y = (fxp1u1 + fyp1u4 + fzp1u7) - Y_factor;
  z = (fxp1u2 + fyp1u5 + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] =  1; ixyz[2] =  1; }

  //if (closest2 != 0.0 && min < closest2) return (min);
//  this->ClosestImage(a1, a2, ixyz);
//  fprintf(stdout,"DEBUG: Predict  = %2i %2i %2i\n",ixyz[0],ixyz[1],ixyz[2]);

//  ix = ixyz[0];
//  iy = ixyz[1];
//  iz = ixyz[2];

//D = sqrt(min);
//  fprintf(stdout,"DEBUG: MinDist  = %2i %2i %2i = %8.3f\n", ixmin, iymin, izmin, D);
//  printf("---------------------------------------------------------------\n");
  return(min);
}

// Frame::DIST2_ImageOrtho()
/** Return the minimum orthorhombic imaged distance^2 between coordinates a1 
  * and a2.
  */
double DIST2_ImageOrtho(double *a1, double *a2, double *box) {
  // If box lengths are zero no imaging possible
  if (box[0]==0.0 && box[1]==0.0 && box[2]==0.0) return -1.0;
  double x = a1[0] - a2[0];
  double y = a1[1] - a2[1];
  double z = a1[2] - a2[2];
  // Get rid of sign info
  if (x<0) x=-x;
  if (y<0) y=-y;
  if (z<0) z=-z;
  // Get rid of multiples of box lengths 
  while (x > box[0]) x = x - box[0];
  while (y > box[1]) y = y - box[1];
  while (z > box[2]) z = z - box[2];
  // Find shortest distance in periodic reference
  double D = box[0] - x;
  if (D < x) x = D;
  D = box[1] - y;
  if (D < y) y = D;  
  D = box[2] - z;
  if (D < z) z = D;

  return (x*x + y*y + z*z);
}

// Frame::DIST2_NoImage()
/** Return distance^2 between coordinates in a1 and a2.
  */
double DIST2_NoImage(double *a1, double *a2) {
  double x,y,z,D;

  x = a1[0] - a2[0];
  y = a1[1] - a2[1];
  z = a1[2] - a2[2];

  x=x*x;
  y=y*y;
  z=z*z;

  //D=sqrt(x + y + z);
  D = x + y + z;

  //fprintf(stdout,"Mask1=%8.3lf %8.3lf %8.3lf Mask2=%8.3lf %8.3lf %8.3lf D=%8.3lf\n",
  //        a1[0],a1[1],a1[2],a2[0],a2[1],a2[2],D);

  return D;
}

// boxToRecip()
/// C-accessible boxToRecip routine for use in ptraj_action.c
extern "C" double boxToRecip(double *box, double *ucell, double *recip) {
  double u12x,u12y,u12z;
  double u23x,u23y,u23z;
  double u31x,u31y,u31z;
  double volume,onevolume;

  ucell[0] = box[0]; // ucell(1,1)
  ucell[1] = 0.0;    // ucell(2,1)
  ucell[2] = 0.0;    // ucell(3,1)
  ucell[3] = box[1]*cos(DEGRAD*box[5]); // ucell(1,2)
  ucell[4] = box[1]*sin(DEGRAD*box[5]); // ucell(2,2)
  ucell[5] = 0.0;                       // ucell(3,2)
  ucell[6] = box[2]*cos(DEGRAD*box[4]);                                         // ucell(1,3)
  ucell[7] = (box[1]*box[2]*cos(DEGRAD*box[3]) - ucell[6]*ucell[3]) / ucell[4]; // ucell(2,3)
  ucell[8] = sqrt(box[2]*box[2] - ucell[6]*ucell[6] - ucell[7]*ucell[7]);       // ucell(3,3)

  // Get reciprocal vectors
  u23x = ucell[4]*ucell[8] - ucell[5]*ucell[7];
  u23y = ucell[5]*ucell[6] - ucell[3]*ucell[8];
  u23z = ucell[3]*ucell[7] - ucell[4]*ucell[6];
  u31x = ucell[7]*ucell[2] - ucell[8]*ucell[1];
  u31y = ucell[8]*ucell[0] - ucell[6]*ucell[2];
  u31z = ucell[6]*ucell[1] - ucell[7]*ucell[0];
  u12x = ucell[1]*ucell[5] - ucell[2]*ucell[4];
  u12y = ucell[2]*ucell[3] - ucell[0]*ucell[5];
  u12z = ucell[0]*ucell[4] - ucell[1]*ucell[3];
  volume=ucell[0]*u23x + ucell[1]*u23y + ucell[2]*u23z;
  onevolume = 1.0 / volume;

  recip[0] = u23x*onevolume;
  recip[1] = u23y*onevolume;
  recip[2] = u23z*onevolume;
  recip[3] = u31x*onevolume;
  recip[4] = u31y*onevolume;
  recip[5] = u31z*onevolume;
  recip[6] = u12x*onevolume;
  recip[7] = u12y*onevolume;
  recip[8] = u12z*onevolume;

  return volume;
}

// calculateDistance2
/// C-accessible distance calc routine for functions in ptraj_actions.c
extern "C" double calculateDistance2(int i, int j, double *x, double *y, double *z,
                                     double *box, double *ucell, double *recip,
                                     double closest2, int noimage)
{
  double A0[3];
  double A1[3];
  
  A0[0] = x[i];
  A0[1] = y[i];
  A0[2] = z[i];
  A1[0] = x[j];
  A1[1] = y[j];
  A1[2] = z[j];

  // NO IMAGE
  if (box == NULL || box[0] == 0.0 || noimage > 0)
    return DIST2_NoImage(A0, A1);
  // ORTHORHOMBIC IMAGING  
  if (box[3] == 90.0 && box[4] == 90.0 && box[5] == 90.0)
    return DIST2_ImageOrtho(A0, A1, box);
  // NON-ORTHORHOMBIC
  return DIST2_ImageNonOrtho(A0, A1, ucell, recip);
}

