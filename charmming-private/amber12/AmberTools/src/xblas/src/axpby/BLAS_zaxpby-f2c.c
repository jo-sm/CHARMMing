
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_zaxpby(int n, const void* alpha, const void* x, int incx, 
   const void* beta, void* y, 
   int incy);


extern void FC_FUNC_(blas_zaxpby,BLAS_ZAXPBY)
  ( int* n, const void* alpha, const void* x, int* incx, const void* beta, void* y, int* incy )
{
	 BLAS_zaxpby(* n,  alpha,   x, * incx, 
    beta,  y, 
   * incy);
}
