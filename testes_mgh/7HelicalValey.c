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
  x[0] = -1.0;
  x[1] = 0.0;
  x[2] = 0.0;
}

void evalf(int n, double *x, double *f, int *flag)
{
  int i;

  *flag = 0;
  double tpi, th, r, t;

  tpi = 8.0 * atan(1.0);
  th = copysign(0.25, x[1]);

  if (x[0] > 0.0)
    th = atan(x[1] / x[0]) / tpi;
  if (x[0] < 0.0)
    th = atan(x[1] / x[0]) / tpi + 0.5;

  r = sqrt(pow(x[0], 2) + pow(x[1], 2));
  t = x[2] - 10.0 * th;
  *f = 100.0 * (pow(t, 2) + pow(r - 1.0, 2)) + pow(x[2], 2);
}

void evalg(int n, double *x, double *g, int *flag)
{
  int i;

  *flag = 0;

  double tpi, th, arg, r, t, s1;

  tpi = 8.0 * atan(1.0);
  th = copysign(0.25, x[1]);

  if (x[0] > 0.0)
    th = atan(x[1] / x[0]) / tpi;
  if (x[0] < 0.0)
    th = atan(x[1] / x[0]) / tpi + 0.5;

  arg = pow(x[0], 2) + pow(x[1], 2);
  r = sqrt(arg);
  t = x[2] - 10.0 * th;
  s1 = 10.0 * t / (tpi * arg);

  g[0] = 200.0 * (x[0] - x[0] / r + x[1] * s1);
  g[1] = 200.0 * (x[1] - x[1] / r - x[0] * s1);
  g[2] = 2.0 * (100.0 * t + x[2]);
}

void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
