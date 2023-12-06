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
 **********************a************************************************/

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

   x[0] = -1.2;
   x[1] = 1.0;
}

void evalf(int n, double *x, double *f, int *flag)
{

   int i;
   *flag = 0;
   *f = pow(10.0 * (x[1] - pow(x[0], 2)), 2) + pow(1.0 - x[0], 2);
}

/***********************************************************************
 **********************************************************************/

void evalg(int n, double *x, double *g, int *flag)
{
   int i;

   *flag = 0;

   double t1 = 10.0 * (x[1] - pow(x[0], 2));
   double t2 = 1.0 - x[0];

   g[0] = 2.0 * (-20.0 * x[0] * t1 - t2);
   g[1] = 20.0 * t1;
}

void proj(int n, double *x, int *flag)
{
   int i;

   *flag = 0;
}