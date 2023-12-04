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

int global_m = 20;
// extern void trace_on(int);

/***********************************************************************
 **********************************************************************/

void inidim(int *n)
{
  /* Set problem data */
  *n = 4;
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
  int j;

  x[0] = 25.0;
  x[1] = 5.0;
  x[2] = -5.0;
  x[3] = -1.0;
}

void evalf(int n, double *x, double *f, int *flag)
{
  *flag = 0;
  *f = 0.0;

  trace_on(1);

  adouble *ax = new adouble[n];
  adouble af = 0.0;
  for (int i = 0; i < n; i++)
    ax[i] <<= x[i];

  for (int i = 1; i <= global_m; ++i)
  {
    adouble d1 = (adouble)i / 5.0;
    adouble d2 = sin(d1);
    adouble t1 = ax[0] + d1 * ax[1] - exp(d1);
    adouble t2 = ax[2] + d2 * ax[3] - cos(d1);
    adouble t = pow(t1, 2) + pow(t2, 2);
    af += pow(t, 2);
  }

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
