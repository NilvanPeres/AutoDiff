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
#include <cfloat>
#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

extern double *l, *u;
int global_m = 15;
double y8[15] = {0.14, 0.18, 0.22, 0.25, 0.29, 0.32, 0.35, 0.39, 0.37, 0.58, 0.73, 0.96, 1.34, 2.10, 4.39};
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
  x[0] = 1.0;
  x[1] = 1.0;
  x[2] = 1.0;
}

void evalf(int n, double *x, double *f, int *flag)
{
  *flag = 0;
  *f = 0.0;

  for (int i = 1; i <= 15; i++)
  {
    double t1 = (double)i;
    double t2 = 16.0 - (double)i;
    double t3 = fmin(t1, t2);

    double d1 = t2 * x[1] + t3 * x[2];

    if (d1 != 0.0)
    {
      *f += pow(y8[i - 1] - (x[0] + (t1 / d1)), 2);
    }
    else
    {
      *f = DBL_MAX;
      *flag = -3;
      return;
    }
  }
}

void evalg(int n, double *x, double *g, int *flag)
{
  *flag = 0;
  g[0] = g[1] = g[2] = 0.0;

  for (int i = 1; i <= 15; i++)
  {
    double t1 = (double)i;
    double t2 = 16.0 - (double)i;
    double t3 = fmin(t1, t2);

    double d1 = t2 * x[1] + t3 * x[2];
    double d2 = pow(d1, 2);

    if (d1 != 0.0)
    {
      double s1 = y8[i - 1] - (x[0] + (t1 / d1));
      double s2 = t1 / d2;

      g[0] -= 2.0 * s1;
      g[1] += 2.0 * s1 * s2 * t2;
      g[2] += 2.0 * s1 * s2 * t3;
    }
    else
    {
      g[0] = g[1] = g[2] = DBL_MAX;
      *flag = -3;
      return;
    }
  }
}

void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
