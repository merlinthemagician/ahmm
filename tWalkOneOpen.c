/*
 * tWalkOneOpen.c
 *
 * MCMC fitting of a Markov model
 * 
 *
 * Ivo Siekmann, 29/11/2010
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef HPC
#include <gsl/gsl_matrix.h>
#else
#include "gsl/gsl_matrix.h"
#endif

#include "matrixArith.h"
#include "matrixExp.h"
#include "matrixIO.h"
#include "mctools.h"
#include "mp_parameter.h"
#include "mp_proposal.h"
#include "mp_mcmc.h"
#include "mp_HMMlikelihoods.h"
#include "mp_HMMutils.h"
#include "likelihood.h"

#define EPS 1e-12

#define OUT stdout
#define ERR stderr

#define NDATA 10000000

/* Iterations */
static int nIter = 1e4;

/* sampling interval */
static double tau=0.05;

/* count of open events, "histogram" of closed events, maximal
   consecutive closed observations */
static int NO, *NC, maxClosed=NDATA, nData=NDATA;

/*estimated closed probability */
static double pC=0;

/* Markov model and proposal */
static gsl_matrix_view q1View, q2View;
static gsl_matrix *Qp1, *Qp2;

/* Returns likelihood for Markov model, given a sequence of open and closed states which is encoded in NO, NC*/
double dPosterior(const parameters *p, int nP,
		  const intparameters *ip, int nIP) {
  int which=getTWalkSampled();
  gsl_matrix *Q=which?&q1View.matrix:&q2View.matrix;
  gsl_matrix *Qp=which?Qp1:Qp2;
  double out=logprNOCdiff1(NO, NC, Q, Qp, tau, maxClosed);

/*   fprintf(OUT, "dPosterior(): Parameter:\n"); */
/*   printParameters(OUT, p, nP+qView.matrix.size1); */
/*   fprintf(OUT, "dPosterior(): Proposal:\n"); */
/*   printProposal(OUT, p, nP+qView.matrix.size1); */

/*   exit(1); */

  return  out;
}

/* Returns likelihood for Markov model, given a sequence of open and closed states which is encoded in NO, NC*/
double dPosteriorNew(const parameters *p, int nP,
		     likelihood *L) {
  int which=getTWalkSampled();
  gsl_matrix *Q=which?&q1View.matrix:&q2View.matrix;
  gsl_matrix *Qp=which?Qp1:Qp2;
  /* double out; */

  return  logprNOC1(NO, NC, Q, Qp, tau, maxClosed, L);
  /* return  out; */
}


/* File for output of stationary distribution */
FILE *statfp;

void outRatesStatDist(FILE *fp, int iteration,
		      const parameters *rates, int nP) {
  int which=getTWalkSampled();
  gsl_matrix *Q=which?&q1View.matrix:&q2View.matrix;

  QandPhiout(fp, iteration, rates, nP, statfp, Q);
}

static char *datafn="NOFILE";

const double allAlpha=1;


int main(int argc, char **argv) {
  int i;
  gsl_rng *r=gsl_rng_alloc(gsl_rng_taus);
  int nStates = 0, seed = 42;
  
  int *qIndices;
  int fullTwalk=1, walkTraverse=0;

  char *modelfile="model/rates.txt";

  double delta =0.1, thresh = 15/* -20 */, *model = NULL, *otherModel=NULL;
  parameters *r1, *r2;
  int nR;
  int lines=60000;

  char *ratesFn="rates.dat", *likeFn="likelihood.dat";
  char *statFn="statDist.dat";
  char *restartFn=NULL;
  FILE *ratesfp, *likefp;


  mp_getArgsOneOpen(argv, argc,
		    &datafn,
		    &modelfile,
		    &nStates,
		    &nIter,
		    &delta,
		    &seed,
		    &restartFn,
		    &thresh,
		    &ratesFn,
		    &statFn,
		    &likeFn);

  NC=columnsData2NOC(datafn, thresh, &NO, &maxClosed,
		     &lines, &nData, &pC);
  /* exit(1); */
  /* Read model from file modelfile in double vector model */
  model=readModel(modelfile, nStates);
  otherModel=malloc(nStates*nStates*sizeof(double));
  memcpy(otherModel, model, nStates*nStates*sizeof(double));

  /* exit(1); */

  setOutput(outRatesStatDist);

  ratesfp=fopen(ratesFn, "w"), likefp=fopen(likeFn, "w");

  statfp=fopen(statFn, "w");
  fprintf(statfp, "Iteration\t");
  for(i=1; i<=nStates; i++) {
    fprintf(statfp, "\"{/Symbol p}_%i\"\t", i);
  }
  fprintf(statfp, "\n");
  fflush(statfp);

  fprintf(OUT, "Initialising rate constants...\n");
  q1View=gsl_matrix_view_array(model, nStates, nStates);
  makeConservative(&q1View.matrix);
  Qp1=gsl_matrix_alloc(nStates, nStates);
  gsl_matrix_memcpy(Qp1, &q1View.matrix);
  nR=countPositive(&q1View.matrix);

  q2View=gsl_matrix_view_array(otherModel, nStates, nStates);
  makeConservative(&q2View.matrix);
  Qp2=gsl_matrix_alloc(nStates, nStates);
  gsl_matrix_memcpy(Qp2, &q2View.matrix);


  /* Parameters: Rates + stationary probabilities */
  r1=malloc(nR*sizeof(parameters));
  r2=malloc(nR*sizeof(parameters));

  fprintf(OUT, "Rates below diagonal first...\n");
  /* Initialise parameters: rates + stationary probabilities */
  initialiseParameters(r1, 0, 0, 10000, nR);
  initialiseParameters(r2, 0, 0, 10000, nR);
  qIndices=calloc(nR, sizeof(int));
  gsl_rng_set(r,seed);
  
  matrix2ParameterSubdiag(&q1View.matrix, Qp1,
			  r1, nR,qIndices);
  matrix2ParameterSubdiag(&q2View.matrix, Qp2,
			  r2, nR,qIndices);
  /* matrix2ParameterPairs(&q1View.matrix, Qp1,  */
  /* 			r1, nR,qIndices); */
  /* matrix2ParameterPairs(&q2View.matrix, Qp2,  */
  /* 			r2, nR,qIndices); */
  fprintf(OUT, "Rates saved to parameters...\n");
  printParameters(OUT,r1,nR);
  printParameters(OUT,r2,nR);

    fprintf(OUT, "Initialising first sample...\n");
    initUniform(r, r1, nR);
    fprintf(OUT, "Initialising second sample...\n");
    for (i=0; i<nR; i++) {
      double rOrg=getParameter(r1, i);
      double rPush=walkPositive(r, rOrg, delta);
      /* fprintf(OUT, "Setting parameter r2[%i] from %f to %f\n", i, rOrg, rOrg+rPush); */
      setProposal(r1, i, rOrg);
      setParameter(r2, i, rOrg+rPush);
      setProposal(r2, i, rOrg+rPush);
    }
    /* printParameters(OUT, r2, nR); */

    fprintf(OUT, "%s", "Done...\n");

    initpMove(nR);

  if(fullTwalk){
    iterateDoubleTwalk(ratesfp, likefp,
		       sampleTwalkPositive,
		       expDoublePrior, dPosteriorNew,
		       r1, nR, r2, nR,
		       seed, nIter);
  }
  else {
    if(walkTraverse)
      iterateDoubleTwalk(ratesfp, likefp,
			 sampleWalkTraversePositivePart,
			 expDoublePrior, dPosteriorNew,
			 r1, nR, r2, nR,
			 seed, nIter);
    else
      iterateDoubleTwalk(ratesfp, likefp,
			 sampleWalkTraversePositivePart,
			 expDoublePrior, dPosteriorNew,
			 r1, nR, r2, nR,
			 seed, nIter);
  }

  return 0;
}
