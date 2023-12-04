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
    x[j] = j + 1;
  }
}

void evalf(int n, double *x, double *f, int *flag)
{
  *flag = 0;
  double t1 = -0.25;
  double t2 = 0.0;

  for (int j = 0; j < n; j++)
  {
    t1 += pow(x[j], 2);
    t2 += pow(x[j] - 1.0, 2);
  }

  *f = 1.0e-05 * t2 + pow(t1, 2);
}

void evalg(int n, double *x, double *g, int *flag)
{
  *flag = 0;

  double s1 = 0.0;
  for (int i = 0; i < n; i++)
  {
    s1 += (i + 1) * x[i];
  }
  double w3[n];
  for (int i = 0; i < n; i++)
  {
    w3[i] = (i + 1) * ((i + 1) * s1 - 1.0);
  }
  double d1 = 0.0;
  for (int i = 0; i < n; i++)
  {
    d1 += w3[i];
  }
  d1 *= 2.0;
  for (int i = 0; i < n; i++)
  {
    g[i] = d1 * (i + 1);
  }
}

void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
