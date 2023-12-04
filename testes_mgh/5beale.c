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
double y5[3] = {1.5, 2.25, 2.625};

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
  *f = 0.0;
  double s1, t1;

  for (int i = 1; i <= 3; i++)
  {
    s1 = 1.0 - pow(x[1], i);
    t1 = y5[i - 1] - x[0] * s1;
    *f += pow(t1, 2);
  }
}

void evalg(int n, double *x, double *g, int *flag)
{
  int i;

  *flag = 0;

  double t1;

  for (int j = 1; j <= 3; j++)
  {
    g[j - 1] = 0.0;
  }

  for (int i = 1; i <= 3; i++)
  {
    t1 = y5[i - 1] - x[0] * (1.0 - pow(x[1], i));

    g[0] = g[0] + 2 * t1 * (pow(x[1], i) - 1.0);
    g[1] = g[1] + 2 * t1 * (i * x[0] * pow(x[1], (i - 1)));
    // g[0] = 2.0 * t1 * (pow(x[1], i) - 1.0);
    // g[1] = 2.0 * t1 * (i * x[0] * pow(x[1], (i - 1)));
  }
}

void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}