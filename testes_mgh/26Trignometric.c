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
  *n = 10;
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
  int j;

  x[0] = 1.0 / (double)n;
}

void evalf(int n, double *x, double *f, int *flag)
{
  *flag = 0;
  double s1 = 0.0;

  for (int j = 0; j < n; ++j)
  {
    s1 += cos(x[j]);
  }

  *f = 0.0;

  for (int j = 0; j < n; ++j)
  {
    double t = (double)(n + j + 1) - sin(x[j]) - s1 - (double)(j + 1) * cos(x[j]);
    *f += pow(t, 2);
  }
}

void evalg(int n, double *x, double *g, int *flag)
{
  *flag = 0;
  double s1 = 0.0;
  double s2 = 0.0;

  for (int j = 0; j < n; ++j)
  {
    g[j] = cos(x[j]);
    s1 += g[j];
  }

  for (int j = 0; j < n; ++j)
  {
    double th = sin(x[j]);
    double t = (double)(n + j + 1) - th - s1 - (double)(j + 1) * g[j];
    s2 += t;
    g[j] = ((double)(j + 1) * th - g[j]) * t;
  }

  for (int j = 0; j < n; ++j)
  {
    g[j] = 2.0 * (g[j] + sin(x[j]) * s2);
  }
}

void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
