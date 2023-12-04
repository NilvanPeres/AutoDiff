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

  x[0] = -1.0;
}

void evalf(int n, double *x, double *f, int *flag)
{
  *flag = 0;

  for (int i = 0; i < n; i++)
  {
    double s1 = x[i] * (2.0 + 5.0 * pow(x[i], 2)) + 1.0;
    for (int j = max(0, i - 5); j < min(n, i + 1); j++)
    {
      if (j != i)
      {
        s1 -= x[j] * (1.0 + x[j]);
      }
    }
    *f += pow(s1, 2);
  }
}

void evalg(int n, double *x, double *g, int *flag)
{
  *flag = 0;

  for (int i = 0; i < n; i++)
  {
    g[i] = 0.0;
  }
  for (int i = 0; i < n; i++)
  {
    double s1 = 0.0;
    for (int j = max(0, i - 5); j < i; j++)
    {
      s1 += x[j] * (1.0 + x[j]);
    }
    if (i != n - 1)
    {
      s1 += x[i + 1] * (1.0 + x[i + 1]);
    }
    double t = x[i] * (2.0 + 5.0 * pow(x[i], 2)) + 1.0 - s1;
    for (int j = max(0, i - 5); j < i; j++)
    {
      g[j] -= 2.0 * t * (1.0 + 2.0 * x[j]);
    }
    g[i] += 2.0 * t * (2.0 + 15.0 * pow(x[i], 2));
    if (i != n - 1)
    {
      g[i + 1] -= 2.0 * t * (1.0 + 2.0 * x[i + 1]);
    }
  }
}
void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
