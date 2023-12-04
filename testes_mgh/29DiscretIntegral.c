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

  double d1 = 1.0 / (n + 1.0);
  double w1[n + 1], w2[n + 1];
  *f = 0.0;
  w1[0] = d1 * pow((x[0] + d1 + 1.0), 3);
  w2[n] = (1.0 - n * d1) * pow((x[n - 1] + n * d1 + 1.0), 3);
  w2[n + 1] = 0.0;
  for (int i = 1; i < n; i++)
  {
    double t1 = i * d1;
    double t2 = (n - i + 1) * d1;
    w1[i] = w1[i - 1] + t1 * pow((x[i] + t1 + 1.0), 3);
    w2[n - i] = w2[n - i + 1] + (1.0 - t2) * pow((x[n - i] + t2 + 1.0), 3);
  }
  for (int i = 0; i < n; i++)
  {
    double t1 = i * d1;
    *f += pow((x[i] + 0.5 * d1 * ((1.0 - t1) * w1[i] + t1 * w2[i + 1])), 2);
  }
}

void evalg(int n, double *x, double *g, int *flag)
{
  *flag = 0;

  double d1 = 1.0 / (n + 1.0);
  double w1[n + 1], w2[n + 1];
  for (int i = 0; i < n; i++)
  {
    g[i] = 0.0;
  }
  w1[0] = d1 * pow((x[0] + d1 + 1.0), 3);
  w2[n] = (1.0 - n * d1) * pow((x[n - 1] + n * d1 + 1.0), 3);
  w2[n + 1] = 0.0;
  for (int i = 1; i < n; i++)
  {
    double t1 = i * d1;
    double t2 = (n - i + 1) * d1;
    w1[i] = w1[i - 1] + t1 * pow((x[i] + t1 + 1.0), 3);
    w2[n - i] = w2[n - i + 1] + (1.0 - t2) * pow((x[n - i] + t2 + 1.0), 3);
  }
  for (int i = 0; i < n; i++)
  {
    double t1 = i * d1;
    double t = x[i] + 0.5 * d1 * ((1.0 - t1) * w1[i] + t1 * w2[i + 1]);
    for (int j = 0; j <= i; j++)
    {
      double t2 = j * d1;
      g[j] += 2.0 * t * (1.5 * d1 * (1.0 - t1) * t2 * pow((x[j] + t2 + 1.0), 2));
    }
    g[i] += 2.0 * t;
    for (int j = i + 1; j < n; j++)
    {
      double t2 = j * d1;
      g[j] += 2.0 * t * (1.5 * d1 * (1.0 - t2) * t1 * pow((x[j] + t2 + 1.0), 2));
    }
  }
}
void proj(int n, double *x, int *flag)
{
  int i;

  *flag = 0;
}
