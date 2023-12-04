/* =================================================================
   File: toyprob.c
   =================================================================

   =================================================================
   Module: Spectral Projected Gradient. Problem definition.
   =================================================================

   Last update of any of the component of this module:

   March 14, 2008.

   Users are encouraged to download periodically updated versions of
   this code at the TANGO Project web page:

   www.ime.usp.br/~egbirgin/tango/

   =================================================================
   ================================================================= */

#include "adolc/adolc.h"
#include "adolc/adouble.h"
#include <math.h>
#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

int global_m = 13;
extern double *l, *u;
// extern void trace_on(int);

/***********************************************************************
 **********************************************************************/

void inidim(int *n)
{
  /* Set problem data */
  *n = 6;
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
  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 1.0;
  x[3] = 1.0;
  x[4] = 1.0;
  x[5] = 1.0;
}

void evalf(int n, double *x, double *f, int *flag)
{
  *flag = 0;
  *f = 0.0;

  trace_on(1);
  adouble *ax = new adouble[n];
  for (int i = 0; i < n; i++)
    ax[i] <<= x[i];
  adouble af = 0.0;

  for (int i = 1; i <= global_m; ++i)
  {
    adouble d1 = (adouble)i / 10.0;
    adouble d2 = exp(-d1) - 5.0 * exp(-10.0 * d1) + 3.0 * exp(-4.0 * d1);
    adouble s1 = exp(-d1 * ax[0]);
    adouble s2 = exp(-d1 * ax[1]);
    adouble s3 = exp(-d1 * ax[4]);
    adouble t = ax[2] * s1 - ax[3] * s2 + ax[5] * s3 - d2;

    af += t * t;
  }
  af >>= *f;
  trace_off();
}

/***********************************************************************
 **********************************************************************/

void evalg(int n, double *x, double *g, int *flag)
{
  *flag = 0;

  gradient(1, n, x, g);
}

/***********************************************************************
 **********************************************************************/

void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
