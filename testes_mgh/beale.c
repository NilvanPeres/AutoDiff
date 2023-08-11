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
  *n = 2;
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
  x[1] = 1.0;
}

void evalf(int n, double *x, double *f, int *flag)
{
  int i;

  *flag = 0;

  // double t1 = 1e+4 * x[0] * x[1] - 1.0;
  // double s1 = exp(-x[0]);
  // double s2 = exp(-x[1]);
  // double t2 = s1 + s2 - 1.0001;

  // *f = pow(t1, 2) + pow(t2, 2);
  double t1, t2, t3;
  double objective = 0.0; // Create a separate variable to accumulate the objective function value

  // for (i = 1; i <= 3; i++)
  // {
  //   s1 = 1.0 - pow(x[1], i);
  //   t1 = 1.5 - x[0] * s1;
  //   objective += pow(t1, 2);
  // }

  t1 = pow(1.5 - x + x * y, 2);
  t2 = pow(2.25 - x + x * pow(y, 2), 2);
  t3 = pow(2.625 - x + x * pow(y, 3), 2);

  *f = term1 + term2 + term3; // Assign the accumulated value to the pointer f
}

void evalg(int n, double *x, double *g, int *flag)
{
  int i;

  *flag = 0;

  // double t1 = 1e+4 * x[0] * x[1] - 1.0;
  // double s1 = exp(-x[0]);
  // double s2 = exp(-x[1]);
  // double t2 = s1 + s2 - 1.0001;

  // g[0] = 2.0 * (1.0e+4 * x[1] * t1 - s1 * t2);
  // g[1] = 2.0 * (1.0e+4 * x[0] * t1 - s2 * t2);

  double t1;
  g[0] = g[1] = 0.0;

  for (int i = 1; i <= 3; i++)
  {
    t1 = 1.5 - x[0] * (1.0 - pow(x[1], i));

    g[0] += 2.0 * t1 * (pow(x[1], i) - 1.0);
    g[1] += 2.0 * t1 * (i * x[0] * pow(x[1], (i - 1)));
  }
}

void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
