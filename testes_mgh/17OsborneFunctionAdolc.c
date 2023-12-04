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
int global_m = 33;

double y17[33] = {8.44E-1, 9.08E-1, 9.32E-1, 9.36E-1, 9.25E-1, 9.08E-1, 8.81E-1, 8.5E-1, 8.18E-1, 7.84E-1,
                  7.51E-1, 7.18E-1, 6.85E-1, 6.58E-1, 6.28E-1, 6.03E-1, 5.8E-1, 5.58E-1, 5.38E-1, 5.22E-1,
                  5.06E-1, 4.9E-1, 4.78E-1, 4.67E-1, 4.57E-1, 4.48E-1, 4.38E-1, 4.31E-1, 4.24E-1, 4.2E-1,
                  4.14E-1, 4.11E-1, 4.06E-1};
/***********************************************************************
 **********************************************************************/

void inidim(int *n)
{
  /* Set problem data */
  *n = 5;
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
  x[0] = 0.5;
  x[1] = 1.5;
  x[2] = -1.0;
  x[3] = 0.01;
  x[4] = 0.02;
}

void evalf(int n, double *x, double *f, int *flag)
{
  *flag = 0;
  *f = 0.0;

  trace_on(1);
  adouble *ax = new adouble[n];
  for (int i = 0; i < n; i++)
    ax[i] <<= x[i];
  adouble af = 0.0;

  for (int i = 1; i <= 33; i++)
  {
    double t1 = 10.0 * (i - 1);

    af += pow(y17[i - 1] - (ax[0] + ax[1] * exp(-t1 * ax[3]) + ax[2] * exp(-t1 * ax[4])), 2);
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
