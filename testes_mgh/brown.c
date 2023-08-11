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
  *n = 2;
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
}

void evalf(int n, double *x, double *f, int *flag)
{
  int i;

  *flag = 0;

  // double t1 = 1e+4 * x[0] * x[1] - 1.0;
  // double s1 = exp(-x[0]);
  // double s2 = exp(-x[1]);
  // double t2 = s1 + s2 - 1.0001;

  // *f = pow(t1, 2) + pow(t2, 2);

  double t1 = x[0] - 1.0e+6;
  double t2 = x[1] - 2.0e-6;
  double t3 = x[0] * x[1] - 2.0;
  *f = pow(t1, 2) + pow(t2, 2) + pow(t3, 2);
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

  double t1 = x[0] - 1.0e+6;
  double t2 = x[1] - 2.0e-6;
  double t3 = x[0] * x[1] - 2.0;
  g[0] = 2.0 * (t1 + x[1] * t3);
  g[1] = 2.0 * (t2 + x[0] * t3);
}

void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
