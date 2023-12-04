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
  x[0] = -3.0;
  x[1] = -1.0;
  x[2] = -3.0;
  x[3] = -1.0;
}

void evalf(int n, double *x, double *f, int *flag)
{
  *flag = 0;

  trace_on(1);

  adouble *ax = new adouble[n];
  for (int i = 0; i < n; i++)
    ax[i] <<= x[i];

  adouble s1 = ax[1] - pow(ax[0], 2);
  adouble s2 = 1.0 - ax[0];
  adouble s3 = ax[1] - 1.0;
  adouble t1 = ax[3] - pow(ax[2], 2);
  adouble t2 = 1.0 - ax[2];
  adouble t3 = ax[3] - 1.0;

  adouble af = 100.0 * pow(s1, 2) + pow(s2, 2) + 90.0 * pow(t1, 2) + pow(t2, 2) + 10.0 * pow(s3 + t3, 2) + pow(s3 - t3, 2) / 10.0;

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
