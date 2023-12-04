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
   *flag = 0;

   trace_on(1);

   adouble *ax = new adouble[n];
   for (int i = 0; i < n; i++)
      ax[i] <<= x[i];

   adouble af = pow(10.0 * (ax[1] - pow(ax[0], 2)), 2) + pow(1.0 - ax[0], 2);

   af >>= *f;

   trace_off();
}

void evalg(int n, double *x, double *g, int *flag)
{
   *flag = 0;

   gradient(1, n, x, g);
}

void proj(int n, double *x, int *flag)
{
   int i;

   *flag = 0;
}