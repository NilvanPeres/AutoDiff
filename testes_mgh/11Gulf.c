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

extern double *l, *u;
int global_m = 99;
// extern void trace_on(int);

/***********************************************************************
 **********************************************************************/

void inidim(int *n)
{
  /* Set problem data */
  *n = 3;
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
  // for ( i=0; i<n; i++ )
  //    x[i] = 60.0;
  x[0] = 5.0;
  x[1] = 2.5;
  x[2] = 1.5e-01;
}

void evalf(int n, double *x, double *f, int *flag)
{
  *flag = 0;

  *f = 0.0;
  double d1 = 2.0 / 3.0;
  double arg, r, t1, t2, t;
  int i;

  for (i = 1; i <= global_m; i++)
  {
    arg = (double)i / 1.0e+2;
    r = fabs(pow(-50.0 * log(arg), d1) + 25.0 - x[1]);
    t1 = pow(r, x[2]) / x[0];
    t2 = exp(-t1);
    t = t2 - arg;
    *f += pow(t, 2);
  }
}

/***********************************************************************
 **********************************************************************/

void evalg(int n, double *x, double *g, int *flag)
{

  *flag = 0;

  double d1 = 2.0 / 3.0;
  double arg, r, t1, t2, t, s1;
  int i;

  for (i = 0; i < 3; i++)
  {
    g[i] = 0.0;
  }

  for (i = 1; i <= global_m; i++)
  {
    arg = (double)i / 1.0e+2;
    r = fabs(pow(-50.0 * log(arg), d1) + 25.0 - x[1]);
    t1 = pow(r, x[2]) / x[0];
    t2 = exp(-t1);
    t = t2 - arg;
    s1 = t1 * t2 * t;
    g[0] = g[0] + s1;
    g[1] = g[1] + s1 / r;
    g[2] = g[2] - s1 * log(r);
  }

  g[0] = 2.0 * g[0] / x[0];
  g[1] = 2.0 * x[2] * g[1];
  g[2] = 2.0 * g[2];
}

/***********************************************************************
 **********************************************************************/

void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
