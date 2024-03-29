/* =================================================================
   File: spg.h
   =================================================================

   =================================================================
   Module: Spectral Projected Gradient. Method subroutine.
   =================================================================

   Last update of any of the component of this module:

   March 14, 2008.

   Users are encouraged to download periodically updated versions of
   this code at the TANGO Project web page:

   www.ime.usp.br/~egbirgin/tango/

   =================================================================
   ================================================================= */

void spg(int n, double *x, double epsopt, int maxit, int maxfc, int iprint,
         double *f, double *gpsupn, int *iter, int *fcnt, int *spginfo,
         int *inform);

void ls(int n, double *x, double f, double *g, double *d, double *lastfv,
        int maxfc, int *fcnt, double *fnew, double *xnew, int *lsinfo,
        int *inform);

void sevalf(int n, double *x, double *f, int *inform);

void evalf(int n, double *x, double *f, int *flag);

void sevalg(int n, double *x, double *g, int *inform);

void evalg(int n, double *x, double *g, int *flag);

void sproj(int n, double *x, int *inform);

void proj(int n, double *x, int *flag);

void reperr(int inform);

void inidim(int *n);

void inip(int n, double *x, double *l, double *u);
