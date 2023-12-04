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
  x[0] = 0.5;
  x[1] = -2;
}

void evalf(int n, double *x, double *f, int *flag)
{

  *flag = 0;

  double t1 = -13.0 + x[0] + (((5.0 - x[1]) * x[1] - 2.0) * x[1]);
  double t2 = -29.0 + x[0] + (((x[1] + 1.0) * x[1] - 14.0) * x[1]);

  *f = pow(t1, 2) + pow(t2, 2);
}
/***********************************************************************
 **********************************************************************/

void evalg(int n, double *x, double *g, int *flag)
{
  int i;

  *flag = 0;

  double t1 = -13.0 + x[0] + ((5.0 - x[1]) * x[1] - 2.0) * x[1];
  double t2 = -29.0 + x[0] + ((x[1] + 1.0) * x[1] - 14.0) * x[1];

  g[0] = 2.0 * (t1 + t2);
  g[1] = 2.0 * (((10.0 - 3.0 * x[1]) * x[1] - 2.0) * t1 + ((3.0 * x[1] + 2.0) * x[1] - 14.0) * t2);
}

/***********************************************************************
 **********************************************************************/

void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
