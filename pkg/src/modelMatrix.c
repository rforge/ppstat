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
#include <Rmath.h>
#include <Rinternals.h>

SEXP computePointProcessFilterMatrix(SEXP t, SEXP B, SEXP delta, SEXP s, SEXP zero){

  /*  We assume that t and s are in increasing order 
   *  t:     time points where filter matrix is evaluated
   *  B:     matrix of basis function evaluations
   *  delta: distance between basis function evaluations in B
             (they are assumed equidistant)
   *  s:     observed points
   */

  int i, j, nt, ns, nss, *nB, lookupIndex, entry, col;
  double *xt, *xs, *xB, *xZ, d, w, diff, antip, target;
  SEXP Z, BDIM;

  if(!isMatrix(B)) error("B must be a matrix when using computePointProcessFilterMatrix");
 
  nt = length(t);
  BDIM = getAttrib(B, R_DimSymbol);
  PROTECT(BDIM = coerceVector(BDIM,INTSXP));
  nB = INTEGER(BDIM);
  PROTECT(Z = allocMatrix(REALSXP,nt,nB[1]));
  PROTECT(t = coerceVector(t,REALSXP));
  PROTECT(B = coerceVector(B,REALSXP));
  PROTECT(delta = coerceVector(delta,REALSXP));
  PROTECT(s = coerceVector(s,REALSXP));
  PROTECT(zero = coerceVector(zero,REALSXP));
  xt = REAL(t);
  xB = REAL(B);
  xs = REAL(s);
  xZ = REAL(Z);

  d = REAL(delta)[0];
  antip = d*(REAL(zero)[0]-1);
  w = d*(nB[0]-1);
  
  for(j = 0; j < nB[1]; j++) {
    ns = 0;
    col =  nB[0]*j;
    for(i = 0; i < nt; i++) {
      entry = i + nt*j;
      xZ[entry] = 0;
      target = xt[i] + antip;
      while(ns < length(s)-1 && target > xs[ns+1]) ns++;
      nss = ns;
      diff = target - xs[ns];
      if(diff > 0) {
	while(diff <= w) 
	  { 
	    lookupIndex = floor(diff/d);
	    xZ[entry] += xB[lookupIndex + col];
	    ns--;
	    if(ns < 0) break;
	    diff = target - xs[ns];
	  }
      }
      ns = nss;
    }
  }

  UNPROTECT(7);
  return(Z);
}

SEXP computeContinuousProcessFilterMatrix(SEXP t, SEXP B, SEXP delta, SEXP s, SEXP zero){

  /*  We assume that t and s are in increasing order
   *  t:     time points where filter matrix is evaluated
   *  B:     matrix of basis function evaluations
   *  delta: distance between basis function evaluations in B
             (they are assumed equidistant)
   *  s:     observed values
   */

  int i, j, k, nt, *nB, lookupIndexLeft, lookupIndexRight, entry, col;
  double *xt, *xs, *xB, *xZ, *BB, d, w, xs0, diffLeft, diffRight;
  SEXP Z, BDIM;

  if(!isMatrix(B)) error("B must be a matrix when using computeModelMatrix.");
  if(length(t) != length(s)) error("Lengths of vectors of observed values and observation points differ.");

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

  BB = (double *) R_alloc(nB[0],sizeof(double));

  d = REAL(delta)[0];
  w = d*(nB[0]-1);

  xs0 = xs[0];
  for(i = 0; i < nt; i++) {
    xs[i] -= xs0;
  } 

  //  printf("Inside function\n");
  
  for(j = 0; j < nB[1]; j++) {
    //printf("Inside outer loop. j: %d\n",j);
    col =  nB[0]*j;
    xZ[nt*j] = xs0;
    // Integration of basis functions, NOT NEEDED ANYWAY
    // Observed process is forced 0 at zero above!
    /*    BB[0] = 0;
	  for(i = 1; i < nB[0]; i++) {
	  BB[i] = BB[i-1] + d*xB[i];
	  }
    */
    // Computation of filter basis
    for(i = 1; i < nt; i++) {
      //printf("Inside inner loop. i: %d\n",i);
      entry = i + nt*j; 
      k = i-1;
      diffLeft = xt[i] - xt[k];
      lookupIndexLeft = floor(diffLeft/d);
      xZ[entry] = xB[lookupIndexLeft + col]*xs[k]*diffLeft;
      diffRight = diffLeft;
      k --;
      // printf("k: %d  diffLeft: %4.2f\n",k,diffLeft);
      while(k >= 0 && (diffLeft = xt[i] - xt[k]) <= w) 
	{ 
	    lookupIndexLeft = floor(diffLeft/d);
	    lookupIndexRight = floor(diffRight/d);
	    // This implementation of numerical integration uses a trapzoidal
	    // rule for the filter function but the left value for the process.
	    // The idea is to mimic 'predictability' and not anticipate 
	    // (numerically) a jump if the process is e.g. a counting process. 
	    xZ[entry] += (xB[lookupIndexLeft + col] + xB[lookupIndexRight + col])*xs[k]*(xt[k+1]-xt[k])/2;
	    diffRight = diffLeft;
	    k--;
	  }
      xZ[entry] += xs0;
    }
  }
  

  UNPROTECT(6);
  return(Z);
}

SEXP computeFilterMatrix(SEXP t, SEXP B, SEXP delta, SEXP s, SEXP zero, SEXP type){
  SEXP Z;
  PROTECT(type = coerceVector(type,STRSXP));
  switch(*CHAR(STRING_ELT(type,0)))
      {
      case 'p':
	PROTECT(Z = computePointProcessFilterMatrix(t, B, delta, s, zero));
	break;
      case 'c':
	PROTECT(Z = computeContinuousProcessFilterMatrix(t, B, delta, s, zero));
	break;
      default:
	PROTECT(Z = allocMatrix(REALSXP,1,1));
      }
  UNPROTECT(2);
  return(Z);
}
