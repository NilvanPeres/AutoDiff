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

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))
#include "adolc/adolc.h"
#include "adolc/adouble.h"

extern double *l, *u;

/***********************************************************************
 **********************************************************************/

void evalf(int n, double *x, double *f, int *flag)
{
   int i;

   *flag = 0;

   *f = 0.0;

   for (i = 0; i < n; i++)
      *f += x[i] * x[i];
}

/***********************************************************************
 **********************************************************************/

// void evalg(int n, double *x, double *g, int *flag)
// {
//    int i;
//    trace_on();

//    *flag = 0;

//    for (i = 0; i < n; i++)
//       g[i] = 2.0 * x[i];
// }

void sevalg(int n, double *x, double *g, int *flga)
{
   *flag = 0.0;
   double *xad = new double[n];
   double *yad = new double[1];

   trace_on(1);

   for (int i = 0; i < n; i++)
   {
      xad[i] = x[i];
   }
   adouble ax;
   ax <<= xad[0];

   adouble af = pow(ax, 3) + 2.0 * pow(ax, 2) + sin(2.0 * ax);

   af >>= yad[0];
   trace_off();
   double grad[n]; // Vetor de derivadas

   gradient(1, n, xad, grad); // Calcular derivadas

   for (int i = 0; i < n; i++)
   {
      g[i] = grad[i];
   }

   delete[] xad;
   delete[] yad;
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
