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
  double h = 1.0 / ((double)n + 1.0);
  for (int j = 0; j < n; j++)
  {
    x[j] = (j + 1) * h * ((j + 1) * h - 1.0);
  }
}

void evalf(int n, double *x, double *f, int *flag)
{
  *flag = 0;
  double d1 = 1.0 / (n + 1.0);
  *f = pow((2.0 * x[0] - x[1] + 0.5 * pow(d1, 2) * pow((x[0] + d1 + 1.0), 3)), 2);
  for (int i = 1; i < n - 1; i++)
  {
    double d2 = i * d1;
    *f += pow((2.0 * x[i] - x[i - 1] - x[i + 1] + 0.5 * pow(d1, 2) * pow((x[i] + d2 + 1.0), 3)), 2);
  }
  double d2 = n * d1;
  *f += pow((2.0 * x[n - 1] - x[n - 2] + 0.5 * pow(d1, 2) * pow((x[n - 1] + d2 + 1.0), 3)), 2);
}

void evalg(int n, double *x, double *g, int *flag)
{
  *flag = 0;

  double d1 = 1.0 / (n + 1.0);
  double d2 = d1;
  double t = 2.0 * x[0] - x[1] + pow(d1, 2) * pow((x[0] + d2 + 1.0), 3) / 2.0;
  g[0] = 2.0 * t * (2.0 + 1.5 * pow(d1, 2) * pow((x[0] + d2 + 1.0), 2));
  g[1] = -2.0 * t;
  for (int i = 1; i < n - 1; i++)
  {
    d2 = i * d1;
    t = 2.0 * x[i] - x[i - 1] - x[i + 1] + pow(d1, 2) * pow((x[i] + d2 + 1.0), 3) / 2.0;
    g[i - 1] -= 2.0 * t;
    g[i] += 2.0 * t * (2.0 + 1.5 * pow(d1, 2) * pow((x[i] + d2 + 1.0), 2));
    g[i + 1] = -2.0 * t;
  }
}

void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
