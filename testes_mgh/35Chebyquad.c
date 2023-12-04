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
double sum(double *arr, int size)
{
  double result = 0.0;
  for (int i = 0; i < size; i++)
  {
    result += arr[i];
  }
  return result;
}
// extern void trace_on(int);

/***********************************************************************
 **********************************************************************/

void inidim(int *n)
{
  /* Set problem data */
  *n = 8;
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
  for (int j = 0; j < n; j++)
  {
    x[j] = (double)(j + 1) / ((double)n + 1.0);
  }
}

void evalf(int n, double *x, double *f, int *flag)
{
  *flag = 0;
  double fvec[n];

  for (int i = 0; i < n; i++)
  {
    fvec[i] = 0.0;
  }

  for (int j = 0; j < n; j++)
  {
    double t1 = 1.0;
    double t2 = 2.0 * x[j] - 1.0;
    double t = 2.0 * t2;

    for (int i = 0; i < n; i++)
    {
      fvec[i] += t2;
      double th = t * t2 - t1;
      t1 = t2;
      t2 = th;
    }
  }

  *f = 0.0;
  double d1 = 1.0 / (double)n;

  for (int i = 0; i < n; i++)
  {
    double t = d1 * fvec[i];

    if (i % 2 == 0)
    {
      t += 1.0 / (pow((double)i, 2) - 1.0);
    }

    *f += t * t;
  }
}

void evalg(int n, double *x, double *g, int *flag)
{
  *flag = 0;
  double w1[n];

  for (int i = 0; i < n; i++)
  {
    w1[i] = 0.0;
  }

  for (int j = 0; j < n; j++)
  {
    double t1 = 1.0;
    double t2 = 2.0 * x[j] - 1.0;
    double t = 2.0 * t2;

    for (int i = 0; i < n; i++)
    {
      w1[i] += t2;
      double th = t * t2 - t1;
      t1 = t2;
      t2 = th;
    }
  }

  double d1 = 1.0 / (double)n;

  for (int i = 0; i < n; i++)
  {
    w1[i] = d1 * w1[i];

    if (i % 2 == 0)
    {
      w1[i] += 1.0 / (pow((double)i, 2) - 1.0);
    }
  }

  for (int j = 0; j < n; j++)
  {
    g[j] = 0.0;
    double t1 = 1.0;
    double t2 = 2.0 * x[j] - 1.0;
    double t = 2.0 * t2;
    double s1 = 0.0;
    double s2 = 2.0;

    for (int i = 0; i < n; i++)
    {
      g[j] += w1[i] * s2;
      double th = 4.0 * t2 + t * s2 - s1;
      s1 = s2;
      s2 = th;
      th = t * t2 - t1;
      t1 = t2;
      t2 = th;
    }
  }

  double d2 = 2.0 * d1;

  for (int i = 0; i < n; i++)
  {
    g[i] = d2 * g[i];
  }
}

void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
