/*
 * mp_HMMlikelihoods.c
 *
 * Likelihood for HMMs
 * 
 * 
 *
 * Ivo Siekmann, 2/11/2010
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */

#include <stdio.h>
#include <string.h>
#include <math.h>

#ifndef HPC
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#else
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_cblas.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#endif

#include "mp_parameter.h"
#include "mctools.h"
#include "matrixIO.h"
#include "matrixArith.h"
#include "matrixExp.h"
#include "likelihood.h"

#define OUT stdout
#define ERR stderr

#define EPS 1e-12

/*log-Probability for n consecutive open states*/
double logPopenSeq(const gsl_matrix *expQ, int n) {
  double expQnn = gsl_matrix_get(expQ, expQ->size1-1,expQ->size2-1);

  /*compute expQnn^n efficiently*/
  return n*log(expQnn);
}

/*log-Probability for n consecutive open states*/
double logPopenSeqDiff(const gsl_matrix *expQ, 
		       const gsl_matrix *expQp, int n) {
  double expQnn = gsl_matrix_get(expQ, expQ->size1-1,expQ->size2-1);
  double expQnnP = gsl_matrix_get(expQp, expQp->size1-1,expQp->size2-1);

  /*compute expQnn^n efficiently*/
  return n*log(expQnnP/expQnn);
}

/*Probability for n consecutive open states*/
double PopenSeq(const gsl_matrix *expQ, int n) {
  return exp(logPopenSeq(expQ,n));
}

/* Compute OC^nO sequence from CO^n and the matrix exponential*/
double PocoSeqCon(const gsl_matrix* COn, const gsl_matrix* expQ) {
  double prC;
  gsl_vector_const_view expQlR=gsl_matrix_const_row(expQ,expQ->size1-1);

  gsl_vector_const_view COnlastC = 
    gsl_matrix_const_column(COn, COn->size2-1);

  gsl_blas_ddot (&expQlR.vector,&COnlastC.vector, &prC);
  return prC;
}

/*Generates projection by left-multiplication of projector, i.e. kill
  the last n lines*/
void generateCO(const gsl_matrix *expQ, gsl_matrix *CO) {
  int j;
  gsl_matrix_memcpy(CO, expQ);
  for(j=0; j<expQ->size2; j++)
    gsl_matrix_set(CO, CO->size1-1, j, 0);
}


/* Computes the likelihood difference of expQp and expQp for a given
 * sequence of consecutive closed and open states. */
double logprNOCdiff(int NO, int *NC, 
		    const gsl_matrix *expQ, 
		    const gsl_matrix *expQp, int nC) {
  int i;
  double pC, pCp;
  double logPdiff, logClosed;
  gsl_matrix *CO = gsl_matrix_alloc(expQ->size1, expQ->size2);
  gsl_matrix *COp = gsl_matrix_alloc(expQp->size1, expQp->size2);

  gsl_matrix *COn[2], *COnP[2];
  COn[0]= gsl_matrix_alloc(expQ->size1, expQ->size2);
  COn[1]= gsl_matrix_alloc(expQ->size1, expQ->size2);

  COnP[0]= gsl_matrix_alloc(expQp->size1, expQp->size2);
  COnP[1]= gsl_matrix_alloc(expQp->size1, expQp->size2);

  generateCO(expQ,CO);
  generateCO(expQp,COp);

  gsl_matrix_set_identity(COn[0]);
  gsl_matrix_set_identity(COnP[0]);

  /*open sequences*/
  logPdiff=logPopenSeqDiff(expQ, expQp, NO);

  /* fprintf(OUT, "logprNOC(): logOpen(%i): %g\n",NO, logPdiff); */

  /*closed sequence*/
  for(i=0; i<nC; i++) {
    multMat(COn[i%2], CO, COn[(i+1)%2]);
    multMat(COnP[i%2], COp, COnP[(i+1)%2]);

    if(NC[i] > 0) {
      pC=PocoSeqCon(COn[(i+1)%2],expQ);
      pCp=PocoSeqCon(COnP[(i+1)%2],expQp);

      logClosed = NC[i]*log(pCp/pC);
      /* fprintf(OUT, "logprNOC(): %i*logClosed(%i): %g\n", NC[i], */
      /* 	      i+1, logClosed); */
      logPdiff += logClosed;
    }
/*       logprob += NC[i]*log(PocoSeq(CO,expQ, i+1)); */
  }
  /* fprintf(OUT, "logP=%f\n",logPdiff); */

  gsl_matrix_free(CO), gsl_matrix_free(COp);
  gsl_matrix_free(COn[0]),gsl_matrix_free(COn[1]);
  gsl_matrix_free(COnP[0]),gsl_matrix_free(COnP[1]);

  return logPdiff;
}

/* Computes the likelihoods of expQp and expQp for a given
 * sequence of consecutive closed and open states. */
double logprNOC(int NO, int *NC, 
		const gsl_matrix *expQ, 
		const gsl_matrix *expQp, int nC, likelihood *L) {
  int i;
  double pC, pCp;
  double logPdiff, logClosed;
  gsl_matrix *CO = gsl_matrix_alloc(expQ->size1, expQ->size2);
  gsl_matrix *COp = gsl_matrix_alloc(expQp->size1, expQp->size2);
  static double L1;
  double L2;

  gsl_matrix *COn[2], *COnP[2];
  COn[0]= gsl_matrix_alloc(expQ->size1, expQ->size2);
  COn[1]= gsl_matrix_alloc(expQ->size1, expQ->size2);

  COnP[0]= gsl_matrix_alloc(expQp->size1, expQp->size2);
  COnP[1]= gsl_matrix_alloc(expQp->size1, expQp->size2);

  generateCO(expQ,CO);
  generateCO(expQp,COp);

  gsl_matrix_set_identity(COn[0]);
  gsl_matrix_set_identity(COnP[0]);

  /*open sequences*/
  /* if(parameterUpdated()) { */
    L1=logPopenSeq(expQ, NO);
  /* } */
  L2=logPopenSeq(expQp, NO);
  logPdiff=logPopenSeqDiff(expQ, expQp, NO);

  /* fprintf(OUT, "logprNOC(): logOpen(%i): %g\n",NO, logPdiff); */

  /*closed sequence*/
  for(i=0; i<nC; i++) {
    multMat(COn[i%2], CO, COn[(i+1)%2]);
    multMat(COnP[i%2], COp, COnP[(i+1)%2]);

    if(NC[i] > 0) {
      /* if ( parameterUpdated()) { */
	pC=PocoSeqCon(COn[(i+1)%2],expQ);
	L1+=NC[i]*log(pC);
      /* } */

      pCp=PocoSeqCon(COnP[(i+1)%2],expQp);
      L2+=NC[i]*log(pCp);
      logClosed = NC[i]*log(pCp/pC);
      /* fprintf(OUT, "logprNOC(): %i*logClosed(%i): %g\n", NC[i], */
      /* 	      i+1, logClosed); */
      logPdiff += logClosed;
    }
/*       logprob += NC[i]*log(PocoSeq(CO,expQ, i+1)); */
  }
  /* fprintf(OUT, "logP=%f\n",logPdiff); */

  gsl_matrix_free(CO), gsl_matrix_free(COp);
  gsl_matrix_free(COn[0]),gsl_matrix_free(COn[1]);
  gsl_matrix_free(COnP[0]),gsl_matrix_free(COnP[1]);

   setPosterior(L, L1,L2);
   return logPdiff;
}

/* Computes the likelihood difference of Q and Qp for a given
 * sequence of consecutive closed and open states. */
double logprNOCdiff1(int NO, int *NC, 
		     gsl_matrix *Q, 
		     gsl_matrix *Qp, double tau, int nC) {
  gsl_matrix *expQ, *expQp;

  makeConservative(Q), makeConservative(Qp);
  expQ=gsl_matrix_alloc(Q->size1, Q->size2);
  expQp=gsl_matrix_alloc(Qp->size1, Qp->size2);
  padeExp(Q, tau, expQ),  padeExp(Qp, tau, expQp);

  /* fprintf(OUT, "Q:\n"); */
  /* printMat(OUT, Q); */
  /* fprintf(OUT, "expQ:\n"); */
  /* printMat(OUT, expQ); */

  /* fprintf(OUT, "Qp:\n"); */
  /* printMat(OUT, Qp); */
  /* fprintf(OUT, "expQp:\n"); */
  /* printMat(OUT, expQp); */
  /* exit(1); */

  return logprNOCdiff(NO, NC, expQ, expQp, nC);
}

/* Computes the likelihoods of Q and Qp for a given
 * sequence of consecutive closed and open states. */
double logprNOC1(int NO, int *NC, 
		     gsl_matrix *Q, 
		     gsl_matrix *Qp, double tau, int nC, likelihood *L) {
  gsl_matrix *expQ, *expQp;

  makeConservative(Q), makeConservative(Qp);
  expQ=gsl_matrix_alloc(Q->size1, Q->size2);
  expQp=gsl_matrix_alloc(Qp->size1, Qp->size2);
  padeExp(Q, tau, expQ),  padeExp(Qp, tau, expQp);

  /* fprintf(OUT, "Q:\n"); */
  /* printMat(OUT, Q); */
  /* fprintf(OUT, "expQ:\n"); */
  /* printMat(OUT, expQ); */

  /* fprintf(OUT, "Qp:\n"); */
  /* printMat(OUT, Qp); */
  /* fprintf(OUT, "expQp:\n"); */
  /* printMat(OUT, expQp); */
  /* exit(1); */

  return logprNOC(NO, NC, expQ, expQp, nC, L);
}
