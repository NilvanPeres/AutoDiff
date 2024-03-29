/* =================================================================
   File: spgma.c
   =================================================================

   =================================================================
   Module: Spectral Projected Gradient. Main program.
   =================================================================

   Last update of any of the component of this module:

   March 14, 2008.

   Users are encouraged to download periodically updated versions of
   this code at the TANGO Project web page:

   www.ime.usp.br/~egbirgin/tango/

   =================================================================
   ================================================================= */

#include <stdlib.h>
#include <stdio.h>
#include <adolc/adolc.h>

#include "spg.h"

/* Global variables that describe the convex set (a box in this
   case). They will be used by the projection subroutine. */

/* Main program */

int main()
{

   int fcnt, iprint, i, inform, iter, maxfc, maxit, n, spginfo;
   double epsopt, f, gpsupn, *x, *l, *u;
   FILE *fp;

   inidim(&n);

   /* Get memory */
   x = (double *)malloc(n * sizeof(double));
   l = (double *)malloc(n * sizeof(double));
   u = (double *)malloc(n * sizeof(double));

   inip(n, x, l, u);

   /* Set solver parameters */

   iprint = 1;
   maxit = 50000;
   maxfc = 10 * maxit;
   epsopt = 1.0e-06;

   /* Call SPG */
   spg(n, x, epsopt, maxit, maxfc, iprint, &f, &gpsupn, &iter, &fcnt, &spginfo, &inform);

   /* Save solution */
   fp = fopen("solution.txt", "w");
   for (i = 0; i < n; i++)
   {
      fprintf(fp, "x[%d] = %e\n", i, x[i]);
      printf(" Pontos finais               : x[%d] = %e\n", i, x[i]);
   }
   // printar soluções

   fclose(fp);

   /* Save statistics */
   fp = fopen("spg-tabline.out", "w");
   fprintf(fp, "%d %d %d %e %e %d\n", n, iter, fcnt, f, gpsupn, spginfo);
   fclose(fp);

   /* Free memory */
   free(x);
   free(l);
   free(u);

   return 0;
}
