/* =================================================================
   File: toyprob.c
   =================================================================

   =================================================================
   Module: Spectral Projected Gradient. Problem definition.
   =================================================================

   Last update of any of the component of this module:

   March 14, 2008.

   Users are encouraged to download periodically updated versions of
   this code at the TANGO Project web page:

   www.ime.usp.br/~egbirgin/tango/

   =================================================================
   ================================================================= */

#include "adolc/adolc.h"
#include "adolc/adouble.h"
#include <math.h>

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

extern double *l, *u;

/***********************************************************************
 **********************************************************************/

void evalf(int n, double *x, double *f, int *flag)
{
   int i;

   *flag = 0;

   f[0] = pow(x[0], 3) + 2 * pow(x[0], 2) + sin(2 * x[0]);
}

/***********************************************************************
 **********************************************************************/

void evalg(int n, double *x, double *g, int *flag)
{
   int i;

   *flag = 0;

   g[0] = 3 * pow(x[0], 2) + 4 * x[0] + 2 * cos(2 * x[0]);
}

/***********************************************************************
 **********************************************************************/

void proj(int n, double *x, int *flag)
{
   int i;

   *flag = 0;

   for (i = 0; i < n; i++)
      x[i] = max(l[i], min(x[i], u[i]));
}
