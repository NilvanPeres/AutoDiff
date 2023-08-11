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
  int i;

  *flag = 0;

  // double t1 = 1e+4 * x[0] * x[1] - 1.0;
  // double s1 = exp(-x[0]);
  // double s2 = exp(-x[1]);
  // double t2 = s1 + s2 - 1.0001;

  // *f = pow(t1, 2) + pow(t2, 2);
  double s1, s2, s3, t1, t2, t3;

  s1 = x[1] - pow(x[0], 2);
  s2 = 1.0 - x[0];
  s3 = x[1] - 1.0;
  t1 = x[3] - pow(x[2], 2);
  t2 = 1.0 - x[2];
  t3 = x[3] - 1.0;

  *f = 100.0 * pow(s1, 2) + pow(s2, 2) + 90.0 * pow(t1, 2) + pow(t2, 2) + 10.0 * pow(s3 + t3, 2) + pow(s3 - t3, 2) / 10.0;
}

// void evalg(int n, double *x, double *g, int *flag)
// {
//   int i;

//   *flag = 0;

//   // double t1 = 1e+4 * x[0] * x[1] - 1.0;
//   // double s1 = exp(-x[0]);
//   // double s2 = exp(-x[1]);
//   // double t2 = s1 + s2 - 1.0001;

//   // g[0] = 2.0 * (1.0e+4 * x[1] * t1 - s1 * t2);
//   // g[1] = 2.0 * (1.0e+4 * x[0] * t1 - s2 * t2);

//   double s1, s2, s3, t1, t2, t3;

//   s1 = x[1] - pow(x[0], 2);
//   s2 = 1.0 - x[0];
//   s3 = x[1] - 1.0;
//   t1 = x[3] - pow(x[2], 2);
//   t2 = 1.0 - x[2];
//   t3 = x[3] - 1.0;

//   g[0] = -2.0 * (2.0e+2 * x[0] * s1 + s2);
//   g[1] = 2.0e+2 * s1 + 20.2 * s3 + 19.8 * t3;
//   g[2] = -2.0 * (1.8e+2 * x[2] * t1 + t2);
//   g[3] = 1.8e+2 * t1 + 20.2 * t3 + 19.8 * s3;
// }

void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
