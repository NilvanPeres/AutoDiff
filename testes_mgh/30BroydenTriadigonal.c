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
  *f = pow((3.0 - 2.0 * x[0]) * x[0] - 2.0 * x[1] + 1.0, 2);
  for (int i = 1; i < n - 1; i++)
  {
    *f += pow((3.0 - 2.0 * x[i]) * x[i] - x[i - 1] - 2.0 * x[i + 1] + 1.0, 2);
  }
  *f += pow((3.0 - 2.0 * x[n - 1]) * x[n - 1] - x[n - 2] + 1.0, 2);
}

void evalg(int n, double *x, double *g, int *flag)
{
  *flag = 0;

  double t = (3.0 - 2.0 * x[0]) * x[0] - 2.0 * x[1] + 1.0;
  g[0] = 2.0 * t * (3.0 - 4.0 * x[0]);
  g[1] = -4.0 * t;

  for (int i = 1; i < n - 1; i++)
  {
    t = (3.0 - 2.0 * x[i]) * x[i] - x[i - 1] - 2.0 * x[i + 1] + 1.0;
    g[i - 1] = g[i - 1] - 2.0 * t;
    g[i] = g[i] + 2.0 * t * (3.0 - 4.0 * x[i]);
    g[i + 1] = -4.0 * t;
  }

  // For i = n
  t = (3.0 - 2.0 * x[n - 1]) * x[n - 1] - x[n - 2] + 1.0;
  g[n - 2] = g[n - 2] - 2.0 * t;
  g[n - 1] = g[n - 1] + 2.0 * t * (3.0 - 4.0 * x[n - 1]);
}
void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
