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
  x[0] = 0.0;
  x[1] = 10.0;
  x[2] = 20.0;
}

void evalf(int n, double *x, double *f, int *flag)
{
  *flag = 0;

  *f = 0.0;
  double d1, d2, s1, s2, s3, t;
  int i;

  for (i = 1; i <= 3; i++)
  {
    d1 = (double)i;
    d2 = d1 / 10.0;
    s1 = exp(-d2 * x[0]);
    s2 = exp(-d2 * x[1]);
    s3 = exp(-d2) - exp(-d1);
    t = s1 - s2 - s3 * x[2];
    *f += pow(t, 2);
  }
}

/***********************************************************************
 **********************************************************************/

void evalg(int n, double *x, double *g, int *flag)
{

  *flag = 0;

  double d1, d2, s1, s2, s3, t, th;
  int i;

  for (i = 0; i < 3; i++)
  {
    g[i] = 0.0;
  }

  for (i = 1; i <= 3; i++)
  {
    d1 = (double)i;
    d2 = d1 / 10.0;
    s1 = exp(-d2 * x[0]);
    s2 = exp(-d2 * x[1]);
    s3 = exp(-d2) - exp(-d1);
    t = s1 - s2 - s3 * x[2];
    th = d2 * t;
    g[0] = g[0] - s1 * th;
    g[1] = g[1] + s2 * th;
    g[2] = g[2] - s3 * t;
  }

  for (i = 0; i < 3; i++)
  {
    g[i] = 2.0 * g[i];
  }
}

/***********************************************************************
 **********************************************************************/

void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
