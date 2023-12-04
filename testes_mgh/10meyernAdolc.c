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
#include <cfloat>
#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

extern double *l, *u;
int global_m = 16;
double y10[16] = {34780.0, 28610.0, 23650.0, 19630.0, 16370.0, 13720.0, 11540.0, 9744.0, 8261.0, 7030.0, 6005.0, 5147.0, 4427.0, 3820.0, 3307.0, 2872.0};
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
  x[0] = 0.02;
  x[1] = 4000.0;
  x[2] = 250.0;
}

void evalf(int n, double *x, double *f, int *flag)
{
  *flag = 0;
  *f = 0.0;

  trace_on(1);
  adouble af = 0.0;
  adouble *ax = new adouble[n];
  for (int i = 0; i < n; i++)
    ax[i] <<= x[i];

  for (int i = 1; i <= global_m; i++)
  {
    adouble t1 = 45.0 + 5.0 * i;
    adouble t2 = t1 + ax[2];

    if (t2 != 0.0)
    {
      af += pow(ax[0] * exp(ax[1] / t2) - y10[i - 1], 2);
    }
    else
    {
      af = INFINITY;
      *flag = -3;
      return;
    }
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
