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
int global_m = 29;
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
  double d1, s1, d2, s2, t;

  *f = 0.0;

  for (int i = 1; i <= global_m; ++i)
  {
    d1 = (double)i / global_m;
    s1 = 0.0;
    d2 = 1.0;

    for (int j = 2; j <= n; ++j)
    {
      s1 += (double)(j - 1) * d2 * x[j - 1];
      d2 *= d1;
    }

    s2 = 0.0;
    d2 = 1.0;

    for (int j = 1; j <= n; ++j)
    {
      s2 += d2 * x[j - 1];
      d2 *= d1;
    }

    t = s1 - pow(s2, 2) - 1.0;
    *f += pow(t, 2);
  }

  double t1 = x[1] - pow(x[0], 2) - 1.0;
  *f += pow(x[0], 2) + pow(t1, 2);
}

void evalg(int n, double *x, double *g, int *flag)
{
  *flag = 0;
  double d1, s1, d2, s2, t, s3;

  for (int i = 0; i < n; ++i)
  {
    g[i] = 0.0;
  }

  for (int i = 1; i <= global_m; ++i)
  {
    d1 = (double)i / global_m;
    s1 = 0.0;
    d2 = 1.0;

    for (int j = 2; j <= n; ++j)
    {
      s1 += (double)(j - 1) * d2 * x[j - 1];
      d2 *= d1;
    }

    s2 = 0.0;
    d2 = 1.0;

    for (int j = 1; j <= n; ++j)
    {
      s2 += d2 * x[j - 1];
      d2 *= d1;
    }

    t = s1 - pow(s2, 2) - 1.0;
    s3 = 2.0 * d1 * s2;
    d2 = 2.0 / d1;

    for (int j = 1; j <= n; ++j)
    {
      g[j - 1] += d2 * ((double)(j - 1) - s3) * t;
      d2 *= d1;
    }
  }

  double t1 = x[1] - pow(x[0], 2) - 1.0;
  g[0] += x[0] * (2.0 - 4.0 * t1);
  g[1] += 2.0 * t1;
}

void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
