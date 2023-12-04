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
int global_m = 15;
double y9[15] = {
    0.0009, 0.0044, 0.0175, 0.0540, 0.1295,
    0.2420, 0.3521, 0.3989, 0.3521, 0.2420,
    0.1295, 0.0540, 0.0175, 0.0044, 0.0009};
// extern void trace_on(int);

/***********************************************************************
 **********************************************************************/

void inidim(int *n)
{
  /* Set problem data */
  *n = 3;
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
  x[0] = 0.4;
  x[1] = 1.0;
  x[2] = 0.0;
}

void evalf(int n, double *x, double *f, int *flag)
{
  *flag = 0;
  trace_on(1);

  adouble *ax = new adouble[n];
  for (int i = 0; i < n; i++)
    ax[i] <<= x[i];
  adouble af = 0.0;

  for (int i = 1; i <= global_m; ++i)
  {
    adouble d1 = 0.5 * (adouble)(i - 1);
    adouble d2 = 3.5 - d1 - ax[2];
    adouble r = exp(-0.5 * ax[1] * pow(d2, 2));
    adouble t = ax[0] * r - y9[i - 1];
    af += pow(t, 2);
  }
  af >>= *f;

  trace_off();
}

void evalg(int n, double *x, double *g, int *flag)
{
  int i;

  *flag = 0;
  gradient(1, n, x, g);
}

void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
