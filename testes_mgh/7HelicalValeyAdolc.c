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

adouble my_copysign(adouble x, adouble y)
{
  return y >= 0 ? fabs(x) : -fabs(x);
}

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
  x[0] = -1.0;
  x[1] = 0.0;
  x[2] = 0.0;
}

void evalf(int n, double *x, double *f, int *flag)
{
  *flag = 0;

  trace_on(1);

  adouble *ax = new adouble[n];
  for (int i = 0; i < n; i++)
    ax[i] <<= x[i];

  adouble tpi = 8.0 * atan(1.0);
  adouble th = my_copysign(0.25, ax[1]);

  if (ax[0] > 0.0)
    th = atan(ax[1] / ax[0]) / tpi;
  if (ax[0] < 0.0)
    th = atan(ax[1] / ax[0]) / tpi + 0.5;

  adouble r = sqrt(pow(ax[0], 2) + pow(ax[1], 2));
  adouble t = ax[2] - 10.0 * th;
  adouble af = 100.0 * (pow(t, 2) + pow(r - 1.0, 2)) + pow(ax[2], 2);

  af >>= *f;

  trace_off();
}

void evalg(int n, double *x, double *g, int *flag)
{
  int i;

  *flag = 0;
  gradient(1, n, x, g);
}

void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
