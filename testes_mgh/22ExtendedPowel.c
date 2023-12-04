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
  *n = 12;
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

  for (j = 0; j < n; j += 4)
  {
    x[j] = 3.0;
    x[j + 1] = -1.0;
    x[j + 2] = 0.0;
    x[j + 3] = 1.0;
  }
}

void evalf(int n, double *x, double *f, int *flag)
{
  *flag = 0;
  *f = 0.0;

  for (int j = 0; j < n; j += 4)
  {
    double t = x[j] + 10.0 * x[j + 1];
    double t1 = x[j + 2] - x[j + 3];
    double s1 = 5.0 * t1;
    double t2 = x[j + 1] - 2.0 * x[j + 2];
    double s2 = pow(t2, 3);
    double t3 = x[j] - x[j + 3];
    double s3 = 10.0 * pow(t3, 3);

    *f += pow(t, 2) + s1 * t1 + s2 * t2 + s3 * t3;
  }
}

void evalg(int n, double *x, double *g, int *flag)
{
  *flag = 0;

  for (int j = 0; j < n; j += 4)
  {
    double t = x[j] + 10.0 * x[j + 1];
    double t1 = x[j + 2] - x[j + 3];
    double s1 = 5.0 * t1;
    double t2 = x[j + 1] - 2.0 * x[j + 2];
    double s2 = 4.0 * pow(t2, 3);
    double t3 = x[j] - x[j + 3];
    double s3 = 20.0 * pow(t3, 3);

    g[j] = 2.0 * (t + s3);
    g[j + 1] = 20.0 * t + s2;
    g[j + 2] = 2.0 * (s1 - s2);
    g[j + 3] = -2.0 * (s1 + s3);
  }
}

void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
