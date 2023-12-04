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
int global_m = 15;
double y9[15] = {
    0.0009, 0.0044, 0.0175, 0.0540, 0.1295,
    0.2420, 0.3521, 0.3989, 0.3521, 0.2420,
    0.1295, 0.0540, 0.0175, 0.0044, 0.0009};
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
  x[0] = 0.4;
  x[1] = 1.0;
  x[2] = 0.0;
}

void evalf(int n, double *x, double *f, int *flag)
{
  *flag = 0;
  *f = 0.0;

  for (int i = 1; i <= global_m; ++i)
  {
    double d1 = 0.5 * (double)(i - 1);
    double d2 = 3.5 - d1 - x[2];
    double r = exp(-0.5 * x[1] * pow(d2, 2));
    double t = x[0] * r - y9[i - 1];
    *f += pow(t, 2);
  }
}

void evalg(int n, double *x, double *g, int *flag)
{
  *flag = 0;

  for (int i = 0; i < n; ++i)
  {
    g[i] = 0.0;
  }

  for (int i = 1; i <= global_m; ++i)
  {
    double d1 = 0.5 * (double)(i - 1);
    double d2 = 3.5 - d1 - x[2];
    double arg = -0.5 * x[1] * pow(d2, 2);
    double r = exp(arg);
    double t = x[0] * r - y9[i - 1];
    double s1 = r * t;
    double s2 = d2 * s1;

    g[0] += 2.0 * s1;
    g[1] -= x[0] * d2 * s2;
    g[2] += 2.0 * x[0] * x[1] * s2;
  }
}

void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
