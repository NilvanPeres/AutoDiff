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
  // for ( i=0; i<n; i++ )
  //    x[i] = 60.0;
  x[0] = 3.0;
  x[1] = -1.0;
  x[2] = 0.0;
  x[3] = 1.0;
}

void evalf(int n, double *x, double *f, int *flag)
{
  *flag = 0;

  trace_on(1);

  adouble af = 0.0;
  adouble *ax = new adouble[n];
  for (int i = 0; i < n; i++)
    ax[i] <<= x[i];

  adouble t1 = ax[0] + 10.0 * ax[1];
  adouble t2 = sqrt(5.0) * (ax[2] - ax[3]);
  adouble t3 = pow((ax[1] - 2.0 * ax[2]), 2);
  adouble t4 = pow(sqrt(10.0) * (ax[0] - ax[3]), 2);
  af = pow(t1, 2) + pow(t2, 2) + pow(t3, 2) + pow(t4, 2);

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
