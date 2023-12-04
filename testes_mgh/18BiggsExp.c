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

  for (int i = 1; i <= global_m; ++i)
  {
    double d1 = (double)i / 10.0;
    double d2 = exp(-d1) - 5.0 * exp(-10.0 * d1) + 3.0 * exp(-4.0 * d1);
    double s1 = exp(-d1 * x[0]);
    double s2 = exp(-d1 * x[1]);
    double s3 = exp(-d1 * x[4]);
    double t = x[2] * s1 - x[3] * s2 + x[5] * s3 - d2;

    *f += t * t;
  }
}

/***********************************************************************
 **********************************************************************/

void evalg(int n, double *x, double *g, int *flag)
{
  *flag = 0;

  for (int i = 0; i < n; i++)
  {
    g[i] = 0.0;
  }

  for (int i = 1; i <= global_m; i++)
  {
    double d1 = (double)i / 10.0;
    double d2 = exp(-d1) - 5.0 * exp(-10.0 * d1) + 3.0 * exp(-4.0 * d1);
    double s1 = exp(-d1 * x[0]);
    double s2 = exp(-d1 * x[1]);
    double s3 = exp(-d1 * x[4]);
    double t = x[2] * s1 - x[3] * s2 + x[5] * s3 - d2;
    double th = d1 * t;

    g[0] -= s1 * th;
    g[1] += s2 * th;
    g[2] += s1 * t;
    g[3] -= s2 * t;
    g[4] -= s3 * th;
    g[5] += s3 * t;
  }

  g[0] *= 2.0 * x[2];
  g[1] *= 2.0 * x[3];
  g[2] *= 2.0;
  g[3] *= 2.0;
  g[4] *= 2.0 * x[5];
  g[5] *= 2.0;
}

/***********************************************************************
 **********************************************************************/

void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
