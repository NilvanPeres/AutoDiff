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
int global_n = 10;
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

  for (j = 0; j < global_n; j += 2)
  {
    x[j] = -1.2;
    x[j + 1] = 1.0;
  }
}

void evalf(int n, double *x, double *f, int *flag)
{
  *flag = 0;
  *f = 0.0;

  for (int j = 0; j < n; j += 2)
  {
    double t1 = 1.0 - x[j];
    double t2 = 10.0 * (x[j + 1] - pow(x[j], 2));
    *f += pow(t1, 2) + pow(t2, 2);
  }
}

void evalg(int n, double *x, double *g, int *flag)
{
  *flag = 0;

  for (int j = 0; j < n; j += 2)
  {
    double t1 = 1.0 - x[j];
    g[j + 1] = 2.0e+2 * (x[j + 1] - pow(x[j], 2));
    g[j] = -2.0 * (x[j] * g[j + 1] + t1);
  }
}

void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
