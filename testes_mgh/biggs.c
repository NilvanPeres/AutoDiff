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
int global_m = 20;
double y5[3] = {1.5, 2.25, 2.625};

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
  int i;

  *flag = 0;

  double res = 0.0;
  double d1, d2, s1, s2, s3, t;

  for (int i = 1; i <= global_m; i++)
  {
    d1 = (double)i / 10.0;
    d2 = exp(-d1) - 5.0 * exp(-10.0 * d1) + 3.0 * exp(-4.0 * d1);
    s1 = exp(-d1 * x[0]);
    s2 = exp(-d1 * x[1]);
    s3 = exp(-d1 * x[4]);
    t = x[2] * s1 - x[3] * s2 + x[5] * s3 - d2;
    res = res + pow(t, 2);
  }

  *f = res;
}

void evalg(int n, double *x, double *g, int *flag)
{
  int i;

  *flag = 0;

  // double t1 = 1e+4 * x[0] * x[1] - 1.0;
  // double s1 = exp(-x[0]);
  // double s2 = exp(-x[1]);
  // double t2 = s1 + s2 - 1.0001;

  // g[0] = 2.0 * (1.0e+4 * x[1] * t1 - s1 * t2);
  // g[1] = 2.0 * (1.0e+4 * x[0] * t1 - s2 * t2);
  double d1, d2, s1, s2, s3, t, th;

  for (int i = 1; i <= global_m; i++)
  {
    d1 = (double)i / 10.0;
    d2 = exp(-d1) - 5.0 * exp(-10.0 * d1) + 3.0 * exp(-4.0 * d1);
    s1 = exp(-d1 * x[0]);
    s2 = exp(-d1 * x[1]);
    s3 = exp(-d1 * x[4]);
    t = x[2] * s1 - x[3] * s2 + x[5] * s3 - d2;
    th = d1 * t;

    g[0] -= s1 * th;
    g[1] += s2 * th;
    g[2] += s1 * t;
    g[3] -= s2 * t;
    g[4] -= s3 * th;
    g[5] += s3 * t;
  }

  g[0] = 2.0 * x[2] * g[0];
  g[1] = 2.0 * x[3] * g[1];
  g[2] = 2.0 * g[2];
  g[3] = 2.0 * g[3];
  g[4] = 2.0 * x[5] * g[4];
  g[5] = 2.0 * g[5];
}

void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
