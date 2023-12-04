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

  for (int i = 1; i <= global_m; i++)
  {
    double t1 = 45.0 + 5.0 * i;
    double t2 = t1 + x[2];

    if (t2 != 0.0)
    {
      *f += pow(x[0] * exp(x[1] / t2) - y10[i - 1], 2);
    }
    else
    {
      *f = INFINITY;
      *flag = -3;
      return;
    }
  }
}

void evalg(int n, double *x, double *g, int *flag)
{
  *flag = 0;

  for (int i = 0; i < n; i++)
    g[i] = 0.0;

  for (int i = 1; i <= global_m; i++)
  {
    double t1 = 45.0 + 5.0 * i;
    double t2 = t1 + x[2];

    if (t2 != 0.0)
    {
      double t3 = exp(x[1] / t2);
      double s1 = x[0] * t3 - y10[i - 1];

      g[0] += 2.0 * s1 * t3;
      g[1] += 2.0 * s1 * (t3 * x[0] / t2);
      g[2] -= 2.0 * s1 * (t3 * x[0] * x[1] / (t2 * t2));
    }
    else
    {
      for (int j = 0; j < n; j++)
        g[j] = INFINITY;
      *flag = -3;
      return;
    }
  }
}

void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
