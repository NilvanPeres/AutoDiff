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

  for (j = 0; j < n; j++)
  {
    x[j] = 0.5;
  }
}

void evalf(int n, double *x, double *f, int *flag)
{
  *flag = 0;
  double t1 = -1.0;
  double t2 = 0.0;
  double t3 = 0.0;
  double d1 = exp(0.1);
  double d2 = 1.0;
  double s2 = 0.0;

  for (int j = 0; j < n; ++j)
  {
    t1 += pow((double)(n - j + 1), 2) * pow(x[j], 2);
    double s1 = exp(x[j] / 10.0);

    if (j != 0)
    {
      double s3 = s1 + s2 - d2 * (d1 + 1.0);
      t2 += pow(s3, 2);
      t3 += pow(s1 - 1.0 / d1, 2);
    }

    s2 = s1;
    d2 = d1 * d2;
  }

  *f = 1.0e-05 * (t2 + t3) + pow(t1, 2) + pow(x[0] - 0.2, 2);
}

void evalg(int n, double *x, double *g, int *flag)
{
  *flag = 0;
  double t1 = -1.0;
  double d1 = exp(0.1);
  double d2 = 1.0;
  double s2 = 0.0;

  for (int j = 0; j < n; ++j)
  {
    t1 += pow((double)(n - j + 1), 2) * pow(x[j], 2);
  }

  double th = 4.0 * t1;

  for (int j = 0; j < n; ++j)
  {
    g[j] = pow((double)(n - j + 1), 2) * x[j] * th;
    double s1 = exp(x[j] / 10.0);

    if (j != 0)
    {
      double s3 = s1 + s2 - d2 * (d1 + 1.0);
      g[j] += 1.0e-5 * s1 * (s3 + s1 - 1.0 / d1) / 5.0;
      g[j - 1] += 1.0e-5 * s2 * s3 / 5.0;
    }

    s2 = s1;
    d2 = d1 * d2;
  }

  g[0] += 2.0 * (x[0] - 0.2);
}

void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
