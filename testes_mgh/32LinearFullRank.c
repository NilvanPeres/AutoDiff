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
  x[0] = 1.0;
}

void evalf(int n, double *x, double *f, int *flag)
{
  *flag = 0;

  double s1 = 0.0;

  for (int i = 0; i < n; i++)
  {
    s1 += x[i];
  }

  double d1 = 2.0 / (double)n;
  *f = 0.0;

  for (int i = 0; i < n; i++)
  {
    double t = x[i] - d1 * s1 - 1.0;
    *f += t * t;
  }

  for (int i = n; i < n; i++)
  {
    double t = -d1 * s1 - 1.0;
    *f += t * t;
  }
}

void evalg(int n, double *x, double *g, int *flag)
{
  *flag = 0;

  double arg = 2.0 / (double)n;
  double s1 = 0.0;

  for (int i = 0; i < n; i++)
  {
    s1 += x[i];
  }

  double w3[n];

  for (int i = 0; i < n; i++)
  {
    w3[i] = x[i] - arg * s1 - 1.0;
  }

  for (int i = n; i < n; i++)
  {
    w3[i] = -arg * s1 - 1.0;
  }

  double s2 = 0.0;

  for (int i = 0; i < n; i++)
  {
    s2 += w3[i];
  }

  for (int i = 0; i < n; i++)
  {
    g[i] = 2.0 * (w3[i] - arg * s2);
  }
}
void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
