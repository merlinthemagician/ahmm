/*
 *  mp_OneOpen.c
 *
 * MCMC fitting of a Markov model
 * 
 *
 * Ivo Siekmann, 10/12/2014
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

static int nStates;

/* sampling interval */
static double tau=0.05;

/* count of open events, "histogram" of closed events, maximal
   consecutive closed observations */
static int NO, *NC, maxClosed=NDATA, nData=NDATA;

/*estimated closed probability */
static double pC=0;


/* Markov model and proposal */
static gsl_matrix_view qView;
static gsl_matrix *Qp;

/* Returns likelihood for Markov model, given a sequence of open and closed states which is encoded in NO, NC*/
double dPosterior(const parameters *p, int nP,
		  const intparameters *ip, int nIP) {
  return  logprNOCdiff1(NO, NC, &qView.matrix, Qp, tau, maxClosed);
}

double dPosteriorNew(const parameters *p, int nP,
		     likelihood *L) {
  return  logprNOC1(NO, NC, &qView.matrix, Qp, tau, maxClosed, L);
}

/* File for output of stationary distribution */
FILE *ratesfp,  *statfp/* , *gibbsfp */;

void outRatesStatDist(FILE *fp, int iteration,
		      const parameters *rates, int nP) {

  QandPhiout(fp, iteration, rates, nP, statfp, &qView.matrix);
}

static char *ratesFn="rates.dat", *statFn="statDist.dat", *likeFn="likelihood.dat";

static char* restartFn=NULL;

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
  int lines=60000;
  int data_nRows=1000, data_nCol=2;
  FILE *datafp;

  /* getArgs(argv, argc); */
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

  /* Read data file */
  /* columnsData2NOC(datafn); */
  datafp=fopen(datafn, "r");
  NC=mp_matrix2OClength(datafp, &data_nRows, data_nCol,
			&NO, &maxClosed,
			0.35, &pC);

  NC=columnsData2NOC(datafn, thresh, &NO, &maxClosed,
		     &lines, &nData, &pC);
  /* exit(1); */
  /* Read model from file modelfile in double vector model */
  model=readModel(modelfile, nStates);
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

  iterateDouble(ratesfp, likefp,
		expDoublePrior, dPosteriorNew,
		rates, nRates,
		seed,  nIter);

  return 0;
}
