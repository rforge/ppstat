/*  Routines for computing the model matrix for the ppstat package.  
 *  These are intended for use with R.
 *
 *     Copyright (C) 2010 Niels Richard Hansen.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 *
 * These functions are distributed in the hope that they will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the
 * GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 *
 */

#include <R.h>
#include <Rinternals.h>

SEXP computeModelMatrix(SEXP t, SEXP B, SEXP delta, SEXP s){

  /*  We assume that t and s are sorted in increasing order 
   *  t:     time points where model matrix is evaluated
   *  B:     matrix of basis function evaluations
   *  delta: distance between basis function evaluations in B
             (they are assumed equidistant)
   *  s:     observed points
   */

  int i, j, nt, ns, nss, *nB, lookupIndex;
  double *xt, *xs, *xB, *xZ, *xdelta, w, diff;
  SEXP Z, BDIM;

  if(!isMatrix(B)) error("B must be a matrix when using computeModelMatrix");
 
  nt = length(t);
  BDIM = getAttrib(B, R_DimSymbol);
  PROTECT(BDIM = coerceVector(BDIM,INTSXP));
  nB = INTEGER(BDIM);
  PROTECT(Z = allocMatrix(REALSXP,nt,nB[1]));
  PROTECT(t = coerceVector(t,REALSXP));
  PROTECT(B = coerceVector(B,REALSXP));
  PROTECT(delta = coerceVector(delta,REALSXP));
  PROTECT(s = coerceVector(s,REALSXP));
  xt = REAL(t);
  xB = REAL(B);
  xs = REAL(s);
  xZ = REAL(Z);

  w = REAL(delta)[0]*nB[0];

  
  for(j = 0; j < nB[1]; j++) {
    ns = 0;
    for(i = 0; i < nt; i++) {
      xZ[i + nt*j] = 0;
      while(ns < length(s)-1 && xt[i] > xs[ns+1]) ns++;
      nss = ns;
      diff = xt[i] - xs[ns];
      // printf("i, t[i], ns, xs[ns], diff: %d %f %d %f %f \n", i, xt[i], ns, xs[ns], diff);
      if(diff > 0) {
	while(diff <= w) 
	  { 
	    lookupIndex = floor(diff/REAL(delta)[0]);
	    xZ[i + nt*j] += xB[lookupIndex + nB[0]*j];
	    ns--;
	    if(ns < 0) break;
	    diff = xt[i] - xs[ns];
	  }
      }
      ns = nss;
    }
  }

  UNPROTECT(6);
  return(Z);
}
