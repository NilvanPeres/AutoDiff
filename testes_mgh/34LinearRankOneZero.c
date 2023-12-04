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
double sum(double *arr, int size)
{
  double result = 0.0;
  for (int i = 0; i < size; i++)
  {
    result += arr[i];
  }
  return result;
}
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

  for (int j = 1; j < n - 1; j++)
  {
    s1 += j + 1 * x[j];
  }

  *f = 2.0;

  for (int i = 1; i < n - 1; i++)
  {
    double t = (i * s1 - 1.0);
    *f += t * t;
  }
}

void evalg(int n, double *x, double *g, int *flag)
{
  *flag = 0;
  double s1 = 0.0;

  for (int i = 2; i < n - 1; i++)
  {
    s1 += i * x[i];
  }

  double w3[n - 2];

  for (int i = 2; i < n - 1; i++)
  {
    w3[i - 2] = (i - 1.0) * ((i - 1.0) * s1 - 1.0);
  }

  double d1 = 2.0 * sum(w3, n - 2);

  g[0] = 0.0;
  g[n - 1] = 0.0;

  for (int i = 2; i < n - 1; i++)
  {
    g[i - 1] = d1 * i;
  }
}
void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
