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
  *n = 40;
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

  x[0] = 0.5;
}

void evalf(int n, double *x, double *f, int *flag)
{
  *flag = 0;
  double s1 = 0.0;
  for (int i = 0; i < n; i++)
  {
    s1 += x[i];
  }
  *f = 0.0;
  for (int i = 0; i < n - 1; i++)
  {
    double t = x[i] + s1 - (n + 1.0);
    *f += pow(t, 2);
  }
  double product = 1.0;
  for (int i = 0; i < n; i++)
  {
    product *= x[i];
  }
  *f += pow((product - 1.0), 2);
}

void evalg(int n, double *x, double *g, int *flag)
{
  *flag = 0;

  double s1 = 0.0;
  double w1[n], w2[n];
  for (int i = 0; i < n; i++)
  {
    g[i] = 0.0;
    s1 += x[i];
  }
  for (int i = 0; i < n - 1; i++)
  {
    double t = x[i] + s1 - (n + 1);
    for (int j = 0; j < n; j++)
    {
      g[j] += 2.0 * t;
    }
    g[i] += 2.0 * t;
  }
  w1[0] = 1.0;
  w2[n - 1] = 1.0;
  for (int i = 1; i < n; i++)
  {
    w1[i] = w1[i - 1] * x[i - 1];
    w2[n - i - 1] = w2[n - i] * x[n - i];
  }
  for (int i = 0; i < n; i++)
  {
    g[i] += 2.0 * (w1[n - 1] * x[n - 1] - 1.0) * w1[i] * w2[i];
  }
}

void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
