/* =================================================================
   File: toyprob.c
   =================================================================

   =================================================================
   Module: Spectral Projected Gradient. Problem definition.
   =================================================================

   =================================================================
   ================================================================= */

#include "adolc/adolc.h"
#include "adolc/adouble.h"
#include <math.h>
#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

extern double *l, *u;
// extern void trace_on(int);

/***********************************************************************
 **********************************************************************/

void inidim(int *n)
{
  /* Set problem data */
  *n = 2;
}

void inip(int n, double *x, double *l, double *u)
{
  /* Define bounds */
  for (int i = 0; i < n; i++)
  {
    l[i] = -1.0e+30;
    u[i] = 1.0e+30;
  }

  /* Define initial Guess */
  x[0] = 0.5;
  x[1] = -2;
}

void evalf(int n, double *x, double *f, int *flag)
{
  int i;

  *flag = 0;

  double t1 = -13.0 + x[0] + (((5.0 - x[1]) * x[1] - 2.0) * x[1]);
  double t2 = -29.0 + x[0] + (((x[1] + 1.0) * x[1] - 14.0) * x[1]);
  *f = pow(t1, 2) + pow(t2, 2);
}

void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
