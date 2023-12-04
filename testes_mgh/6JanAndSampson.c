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
int global_m = 10;
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
  // for ( i=0; i<n; i++ )
  //    x[i] = 60.0;
  x[0] = 0.2578;
  x[1] = 0.2578;
}

void evalf(int n, double *x, double *f, int *flag)
{
  *flag = 0;

  double t1 = exp(x[0]);
  double t2 = exp(x[1]);
  double d1 = t1;
  double d2 = t2;
  *f = 0.0;

  for (int i = 1; i <= global_m; i++)
  {
    *f += pow(2.0 + 2.0 * i - (d1 + d2), 2);
    d1 *= t1;
    d2 *= t2;
  }
}

/***********************************************************************
 **********************************************************************/

void evalg(int n, double *x, double *g, int *flag)
{
  *flag = 0;
  g[0] = 0.0;
  g[1] = 0.0;

  double t1 = exp(x[0]);
  double t2 = exp(x[1]);
  double d1 = t1;
  double d2 = t2;

  for (int i = 1; i <= global_m; i++)
  {
    double s1 = 2.0 * i + 2.0 - (d1 + d2);
    g[0] -= 2.0 * i * s1 * d1;
    g[1] -= 2.0 * i * s1 * d2;
    d1 *= t1;
    d2 *= t2;
  }
}

/***********************************************************************
 **********************************************************************/

void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
