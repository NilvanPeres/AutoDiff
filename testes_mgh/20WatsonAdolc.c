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

  x[0] = 0.0;
}

void evalf(int n, double *x, double *f, int *flag)
{
  *flag = 0;

  trace_on(1);

  adouble *ax = new adouble[n];
  for (int i = 0; i < n; i++)
    ax[i] <<= x[i];

  adouble af = 0.0;
  adouble d1, s1, d2, s2, t;

  for (int i = 1; i <= 29; ++i)
  {
    d1 = (double)i / 29.0;
    s1 = 0.0;
    d2 = 1.0;

    for (int j = 2; j <= n; ++j)
    {
      s1 += (double)(j - 1) * d2 * ax[j - 1];
      d2 *= d1;
    }

    s2 = 0.0;
    d2 = 1.0;

    for (int j = 1; j <= n; ++j)
    {
      s2 += d2 * ax[j - 1];
      d2 *= d1;
    }

    t = s1 - pow(s2, 2) - 1.0;
    af += pow(t, 2);
  }

  adouble t1 = ax[1] - pow(ax[0], 2) - 1.0;
  af += pow(ax[0], 2) + pow(t1, 2);

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
