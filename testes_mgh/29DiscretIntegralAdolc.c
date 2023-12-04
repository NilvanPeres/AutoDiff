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

  trace_on(1);
  adouble d1 = 1.0 / (n + 1.0);
  adouble *w1 = new adouble[n + 2];
  adouble *w2 = new adouble[n + 2];
  adouble *ax = new adouble[n];
  for (int i = 0; i < n; i++)
    ax[i] <<= x[i];
  adouble af = 0.0;
  w1[0] = d1 * pow((ax[0] + d1 + 1.0), 3);
  w2[n] = (1.0 - n * d1) * pow((ax[n - 1] + n * d1 + 1.0), 3);
  w2[n + 1] = 0.0;
  for (int i = 1; i < n; i++)
  {
    adouble t1 = i * d1;
    adouble t2 = (n - i + 1) * d1;
    w1[i] = w1[i - 1] + t1 * pow((ax[i] + t1 + 1.0), 3);
    w2[n - i] = w2[n - i + 1] + (1.0 - t2) * pow((ax[n - i] + t2 + 1.0), 3);
  }
  for (int i = 0; i < n; i++)
  {
    adouble t1 = i * d1;
    af += pow((ax[i] + 0.5 * d1 * ((1.0 - t1) * w1[i] + t1 * w2[i + 1])), 2);
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
