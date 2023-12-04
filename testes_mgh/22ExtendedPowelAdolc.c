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

  trace_on(1);
  adouble *ax = new adouble[n];
  adouble af;
  for (int i = 0; i < n; i++)
    ax[i] <<= x[i];

  for (int j = 0; j < n; j += 4)
  {
    adouble t = ax[j] + 10.0 * ax[j + 1];
    adouble t1 = ax[j + 2] - ax[j + 3];
    adouble s1 = 5.0 * t1;
    adouble t2 = ax[j + 1] - 2.0 * ax[j + 2];
    adouble s2 = pow(t2, 3);
    adouble t3 = ax[j] - ax[j + 3];
    adouble s3 = 10.0 * pow(t3, 3);

    af += pow(t, 2) + s1 * t1 + s2 * t2 + s3 * t3;
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
