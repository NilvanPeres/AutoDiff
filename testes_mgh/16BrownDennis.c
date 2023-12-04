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

int global_m = 20;
// extern void trace_on(int);

/***********************************************************************
 **********************************************************************/

void inidim(int *n)
{
  /* Set problem data */
  *n = 4;
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

  x[0] = 25.0;
  x[1] = 5.0;
  x[2] = -5.0;
  x[3] = -1.0;
}

void evalf(int n, double *x, double *f, int *flag)
{
  *flag = 0;
  *f = 0.0;

  for (int i = 1; i <= global_m; ++i)
  {
    double d1 = (double)i / 5.0;
    double d2 = sin(d1);
    double t1 = x[0] + d1 * x[1] - exp(d1);
    double t2 = x[2] + d2 * x[3] - cos(d1);
    double t = pow(t1, 2) + pow(t2, 2);
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
    double d1 = (double)i / 5.0;
    double s1 = exp(d1);
    double s2 = sin(d1);
    double s3 = cos(d1);
    double t1 = x[0] + d1 * x[1] - s1;
    double t2 = x[2] + x[3] * s2 - s3;
    double t = pow(t1, 2) + pow(t2, 2);

    g[0] += 4.0 * t * (x[0] + d1 * x[1] - s1);
    g[1] += 4.0 * t * (d1 * (x[0] + d1 * x[1] - s1));
    g[2] += 4.0 * t * (x[2] + s2 * x[3] - s3);
    g[3] += 4.0 * t * (s2 * (x[2] + s2 * x[3] - s3));
  }
}

void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
