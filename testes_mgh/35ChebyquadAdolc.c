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
  *n = 8;
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
  for (int j = 0; j < n; j++)
  {
    x[j] = (double)(j + 1) / ((double)n + 1.0);
  }
}

void evalf(int n, double *x, double *f, int *flag)
{
  *flag = 0;

  trace_on(1);
  adouble *ax = new adouble[n];
  adouble af;
  for (int i = 0; i < n; i++)
    ax[i] <<= x[i];

  adouble fvec[n];

  for (int i = 0; i < n; i++)
  {
    fvec[i] = 0.0;
  }

  for (int j = 0; j < n; j++)
  {
    adouble t1 = 1.0;
    adouble t2 = 2.0 * ax[j] - 1.0;
    adouble t = 2.0 * t2;

    for (int i = 0; i < n; i++)
    {
      fvec[i] += t2;
      adouble th = t * t2 - t1;
      t1 = t2;
      t2 = th;
    }
  }

  af = 0.0;
  adouble d1 = 1.0 / (adouble)n;

  for (int i = 0; i < n; i++)
  {
    adouble t = d1 * fvec[i];

    if (i % 2 == 0)
    {
      t += 1.0 / (pow((adouble)i, 2) - 1.0);
    }

    af += t * t;
  }
  af >>= *f;
  trace_off();
}

void evalg(int n, double *x, double *g, int *flag)
{
  *flag = 0;
  gradient(1, n, x, g);
}

void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
