/*
 *  mp_NOpen.c
 *
 * MCMC fitting of a Markov model with arbitrary number of open and
 * closed states.
 * 
 *
 * Ivo Siekmann, 11/02/2015
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
#include <gsl/gsl_randist.h>
#else
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_randist.h"
#endif

#include "matrixArith.h"
#include "matrixExp.h"
#include "matrixIO.h"
#include "mctools.h"
#include "mp_proposal.h"
#include "mp_parameter.h"
#include "mp_mcmc.h"
#include "mp_HMMlikelihoods.h"
#include "mp_HMMutils.h"

#include "nw_data.h"

#define OUT stdout
#define ERR stderr

#define NDATA 10000000

#define MAXCYCLE 10

/* Iterations */
static int nIter = 1e4;

static int nStates, nOpen;

/* sampling interval */
static double tau=0.05;

/* Sequence of open/closed determined by thresholding */
static int *events, nEvents;


/* Markov model and proposal */
static gsl_matrix_view qView;
static gsl_matrix *Qp;

/* Filter for computing state sequence */
static gsl_matrix *filter;

/*inactivates states from the "wrong" class and rescales pS - returns
  scaling factor.*/
double inactivateStates(gsl_vector *pS, int class, int nOpen) {
  int n = pS->size;
  int nClosed = n-nOpen;
  /*if closed (0), deactivate open */
  int i0 = class ? 0: nClosed;
  int nDeact = class ? nClosed: nOpen;
  gsl_vector_view inact = gsl_vector_subvector(pS, i0, nDeact);
  
  gsl_vector_set_zero(&inact.vector);
  return toProbVector(pS);
}



/* next filter step - returns scaling factor (important for
   log-Likelihood computation  */
double forwardFilter(gsl_vector *pS, 
		     const gsl_vector *filter, 
		     const gsl_matrix *expQ, 
		     int c, int cNext, 
		     int nOpen) {
  int n=pS->size;
  int nClosed = n-nOpen;
  int i, j;
  int j0 = c ? nClosed : 0, i0 = cNext ? nClosed : 0;
  int jSteps = c ? nOpen : nClosed, iSteps = cNext ? nOpen : nClosed;

  gsl_vector_set_zero(pS);
  /* THIS EFFICIENT BIT OF CODE CAN'T BE USED ANYMORE - THE SUM IS
     NEEDED FOR LIKELIHOOD COMPUTATION */
/*   if(iSteps == 1) { */
/*     gsl_vector_set(pS, i0, 1); */
/*     return; */
/*   } */
/****************/

/*   fprintf(OUT, "forwardFilter():\n"); */
  for(i=0; i<iSteps; i++) {
    double sum = 0;
    for(j=0; j< jSteps; j++) {
      double pJ = gsl_vector_get(filter, j0+j);
/*       fprintf(OUT, "i0 = %i, i=%i, j0=%i, j=%i, sum = %f\n", */
/* 	      i0, i, j0,j, sum); */
      double expQji = gsl_matrix_get(expQ, j0+j, i0+i);
      sum += expQji * pJ; 
/*       if(!c && !cNext) */
/* 	fprintf(OUT, "i0 = %i, i=%i, j0=%i, j=%i, pJ = %f, expQji = %f, sum = %f\n", */
/* 	      i0, i, j0,j, pJ, expQji, sum); */
    }
/*       if(!c && !cNext) { */
/* /\* 	fprintf(OUT, "i0 = %i, i=%i, j0=%i, j=%i, pI = %f, expQij = %f, sum = %f\n", *\/ */
/* /\* 		i0, i, j0,j, pI, expQij, sum); *\/ */
/* 	fprintf(OUT, "i0 = %i, i=%i, j0=%i, j=%i, sum = %f\n", */
/* 		i0, i, j0,j, sum); */
/* 	printMat(OUT,expQ); */
/*       } */
    gsl_vector_set(pS, i0+i, sum);

/*     if(!c && !cNext) { */
/*       vprint(OUT,pS); */
/*     } */
  }

  return  toProbVector(pS);
/*   if(!c && !cNext) { */
/*     vprint(OUT,pS); */
/*     exit(1); */
/*   } */
}


/*Computing filter --> FORWARD step, computes log-likelihood of total
  sequence */
double classSToFilter(const int *classS, int nClassS, 
		    const gsl_vector *lambda,
		    const gsl_matrix *expQ,
		    int nOpen,
		    gsl_matrix *filter) {
  int n=lambda->size;
  int i;
  gsl_vector *pS =gsl_vector_alloc(n);
  double L;
  double p;

  /*initialising */
/*   gsl_matrix_set_zero(filter); */

  /* fprintf(OUT, "classSToFilter(): First filter element (n=%i)...\n", */
  /* 	  pS->size); */
  /* fprintf(OUT, "filter size: %i, %i\n", filter->size1, filter->size2); */
  gsl_vector_memcpy(pS, lambda);
  p=inactivateStates(pS, classS[0], nOpen);
  L=log(p);
  /* vprint(OUT,pS); */
  /* fprintf(OUT, "L= %f\n", L); */
  
  /* fprintf(OUT, "Computing filter...\n"); */
  gsl_matrix_set_row(filter, 0, pS);
  /* fprintf(OUT, "classSToFilter(): Set first filter element...\n"); */
  /* vprint(OUT, pS); */

  for(i=1; i<nClassS; i++) {
    /* Filter value which was just computed */
    gsl_vector_const_view oldFilterView =  gsl_matrix_const_row(filter, i-1);
    L += log(
	     forwardFilter(pS, &oldFilterView.vector, 
			   expQ, classS[i-1], classS[i], nOpen)
	     );
    /* fprintf(OUT, "L= %f\n", L); */
    gsl_matrix_set_row(filter, i, pS);
  }

  gsl_vector_free(pS);

  return L;

}

/* POSTERIOR */
/* simple Metropolis-Hastings algorithm for n open states */
double logNOpenDiff(const int *events, int nEvents,
		    gsl_matrix *Q, gsl_matrix *Qp, 
		    double tau, int nOpen) {

  gsl_matrix *expQ1=gsl_matrix_alloc(Q->size1, Q->size2);
  gsl_matrix *expQ2=gsl_matrix_alloc(Qp->size1, Qp->size2);
  gsl_vector *phi1 = gsl_vector_alloc(Q->size2);
  gsl_vector *phi2 = gsl_vector_alloc(Qp->size2);
  static double L1;
  double L2;

  makeConservative(Q), makeConservative(Qp);
/*   fprintf(OUT, "Compute stationary distributions...\n"); */
  statMarkovVector(Q, phi1), statMarkovVector(Qp, phi2);
/*   fprintf(OUT, "Compute matrix exponentials...\n"); */

/*   fprintf(OUT, "Compute likelihoods\n"); */
  if(parameterUpdated()) {
    padeExp(Q, tau, expQ1);
    L1 = classSToFilter(events, nEvents, phi1, expQ1, nOpen, filter);
  }

  if(gsl_matrix_max (Qp) < 100) {
    padeExp(Qp, tau, expQ2);
    L2 = classSToFilter(events, nEvents, phi2, expQ2, nOpen, filter);
  }
  else L2 =  -10000;

/*   fprintf(OUT, "logNOpenDiff(): Q:\n"); */
/*   printMat(OUT, Q); */
/*   fprintf(OUT, "phi: "); */
/*   vprint(OUT, phi1); */

/*   fprintf(OUT, "logNOpenDiff(): Qp:\n"); */
/*   printMat(OUT, Qp); */
/*   fprintf(OUT, "phi: "); */
/*   vprint(OUT, phi2); */

/*   fprintf(OUT, "logNOpenDiff(): L1 = %f, L2 = %f, diff = %f\n", L1, L2, L2-L1); */

/*   exit(1); */
  return L2-L1;
}

/* POSTERIOR */
/* simple Metropolis-Hastings algorithm for n open states */
double logNOpen(const int *events, int nEvents,
		gsl_matrix *Q, gsl_matrix *Qp, 
		double tau, int nOpen, likelihood *L) {

  gsl_matrix *expQ1=gsl_matrix_alloc(Q->size1, Q->size2);
  gsl_matrix *expQ2=gsl_matrix_alloc(Qp->size1, Qp->size2);
  gsl_vector *phi1 = gsl_vector_alloc(Q->size2);
  gsl_vector *phi2 = gsl_vector_alloc(Qp->size2);
  static double L1;
  double L2;

  makeConservative(Q), makeConservative(Qp);
/*   fprintf(OUT, "Compute stationary distributions...\n"); */
  statMarkovVector(Q, phi1), statMarkovVector(Qp, phi2);
/*   fprintf(OUT, "Compute matrix exponentials...\n"); */

/*   fprintf(OUT, "Compute likelihoods\n"); */
  if(parameterUpdated()) {
    padeExp(Q, tau, expQ1);
    L1 = classSToFilter(events, nEvents, phi1, expQ1, nOpen, filter);
  }

  if(gsl_matrix_max (Qp) < 100) {
    padeExp(Qp, tau, expQ2);
    L2 = classSToFilter(events, nEvents, phi2, expQ2, nOpen, filter);
  }
  else L2 =  -10000;

  /* fprintf(OUT, "logNOpenDiff(): Q:\n"); */
  /* printMat(OUT, Q); */
  /* fprintf(OUT, "phi: "); */
  /* vprint(OUT, phi1); */

  /* fprintf(OUT, "logNOpenDiff(): Qp:\n"); */
  /* printMat(OUT, Qp); */
  /* fprintf(OUT, "phi: "); */
  /* vprint(OUT, phi2); */

  /* fprintf(OUT, "logNOpenDiff(): L1 = %f, L2 = %f, diff = %f\n", L1, L2, L2-L1); */

/*   exit(1); */
  setPosterior(L, L1,L2);
  return L2-L1;
}

/* Returns likelihood for Markov model, given a sequence of open and closed states which is encoded in NO, NC*/
double dPosterior(const parameters *p, int nP,
		  const intparameters *ip, int nIP) {

  double out=logNOpenDiff(events, nEvents, &qView.matrix, Qp, tau, nOpen);

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

  double out=logNOpen(events, nEvents, &qView.matrix, Qp, tau, nOpen, L);

/*   fprintf(OUT, "dPosterior(): Parameter:\n"); */
/*   printParameters(OUT, p, nP+qView.matrix.size1); */
/*   fprintf(OUT, "dPosterior(): Proposal:\n"); */
/*   printProposal(OUT, p, nP+qView.matrix.size1); */

/*   exit(1); */

  return  out;
}

/* File for output of stationary distribution */
FILE *ratesfp,  *statfp/* , *gibbsfp */;

void outRatesStatDist(FILE *fp, int iteration,
		      const parameters *rates, int nP) {

  QandPhiout(fp, iteration, rates, nP, statfp, &qView.matrix);
}

static char *ratesFn="rates.dat", *statFn="statDist.dat", *likeFn="likelihood.dat";

static char *datafn="NOFILE";

static  char *modelfile="model/rates.txt";

static double delta = 0.1, thresh = 15, *model = NULL;

static int seed = 42;

int main(int argc, char **argv) {
  gsl_rng *r=gsl_rng_alloc(gsl_rng_taus);

  parameters *rates;
  int *qIndices;
  int i;
  int nRates;
  FILE *likefp;
  int nRows=60000;

  /* getArgs(argv, argc); */
  mp_getArgsNOpen(argv, argc,
		  &datafn,
		  &modelfile,
		  &nStates,
		  &nOpen,
		  &nIter,
		  &delta,
		  &seed,
		  &thresh,
		  &ratesFn,
		  &statFn,
		  &likeFn);

  /* Read data file */
  events=mp_columnsData2Events(datafn, thresh, &nRows, &nEvents);
 /* exit(1); */
  /* Read model from file modelfile in double vector model */
  model=readModel(modelfile, nStates);
  filter=gsl_matrix_alloc(nEvents,nStates);
  /* exit(1); */
  setOutput(outRatesStatDist);

  /* NC=calloc(NDATA, sizeof(int)); */
  /* maxClosed=processDataOneOpen(stdin, NC, &NO, &nData, thresh, &pC); */
  
  ratesfp=fopen(ratesFn, "w"), likefp=fopen(likeFn, "w");

  statfp=fopen(statFn, "w");
  fprintf(statfp, "Iteration\t");
  for(i=1; i<=nStates; i++) {
    fprintf(statfp, "\"{/Symbol p}_%i\"\t", i);
  }
  fprintf(statfp, "\n");
  fflush(statfp);


  fprintf(OUT, "Initialising rate constants...\n");
  qView=gsl_matrix_view_array(model, nStates, nStates);
  makeConservative(&qView.matrix);
  Qp=gsl_matrix_alloc(nStates, nStates);
  gsl_matrix_memcpy(Qp, &qView.matrix);
  nRates=countPositive(&qView.matrix);
  qIndices=calloc(nRates, sizeof(int));
  gsl_rng_set(r,seed);

  fprintf(OUT, "Rates below diagonal first...\n");
  rates=malloc(nRates*sizeof(parameters));
  initialiseParameters(rates, 
		       delta, 0, 10000, nRates);
  matrix2ParameterPairs(&qView.matrix, Qp, 
			rates, nRates,qIndices);
  /* printParameters(OUT, rates, nRates); */
  initUniform(r, rates, nRates);
  /* printParameters(OUT, rates, nRates); */
  /* exit(1); */

  iterateDouble(ratesfp, likefp,
		expDoublePrior, dPosteriorNew,
		rates, nRates,
		seed,  nIter);

  return 0;
}
