/*
 * sojournDist.c
 *
 * Computes distribution over an aggregate of states in an aggregated
 * Markov model.
 *
 * Ivo Siekmann, 15/09/16
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>
#include "mctools.h"
#include "matrixIO.h"

#define OUT stdout
#define ERR stderr

/* Dimension of generator, number of states in the aggregate for which
   dwell-time distribution is to be calculated */
static int nStates, nDist;

static gsl_matrix* Q;

/* Print exponential term */
void printExp(double a, double l) {
  fprintf(OUT, "%g * exp(%f * t)", a, l); 
}

/* Print sum of exponential distributions */
void printDist(const gsl_vector *coeff, const gsl_vector *lambda) {
  int i;
  double c=gsl_vector_get(coeff,0),l=gsl_vector_get(lambda, 0);

  printExp(c, l);
  for(i=1; i<coeff->size; i++) {
    c = gsl_vector_get(coeff, i);
    l= gsl_vector_get(lambda, i);
    if(c >= 0)
      fprintf(OUT, " + ");
    printExp(c, l);
  }
  fprintf(OUT, "\n");
}

/* Edges of M are updated according to the vertex permutation p */
void vertexPerm(const gsl_matrix *M, int *p, gsl_matrix *R) {
  int i, j;

  for (i=0; i<M->size1; i++) {
    for(j=0; j<M->size2; j++) {
      int newI=p[i], newJ=p[j];
      double mIJ=gsl_matrix_get(M,i,j);
      gsl_matrix_set(R, newI,newJ,mIJ);
    }
  }
}

/* Converts lambda to time constants */
gsl_vector* timeConst(const gsl_vector * lambda){
  gsl_vector *tc=gsl_vector_alloc(lambda->size);
  int i;

  for(i=0; i<tc->size; i++) {
    double lI=gsl_vector_get(lambda, i);
    gsl_vector_set(tc, i, 1.0/lI);
  }

  return tc;
}

/* Converts lambda to time constants and squares */
gsl_vector* timeConstSqr(const gsl_vector * lambda){
  gsl_vector *tc=timeConst(lambda);
  int i;

  for(i=0; i<tc->size; i++) {
    double lI=gsl_vector_get(tc, i);
    gsl_vector_set(tc, i, lI*lI);
  }

  return tc;
}

/* Converts lambda to time constants and cubes */
gsl_vector* timeConstCube(const gsl_vector * lambda){
  gsl_vector *tc=timeConst(lambda);
  int i;

  for(i=0; i<tc->size; i++) {
    double lI=gsl_vector_get(tc, i);
    gsl_vector_set(tc, i, lI*lI*lI);
  }

  return tc;
}

/* Returns mean of mix of exponential distributions */
double meanExpMix(const gsl_vector* coeff, const gsl_vector* lambda) {
  gsl_vector *v=gsl_vector_alloc(lambda->size);
  gsl_vector *tc=timeConstSqr(lambda);
  gsl_vector_memcpy(v, tc);
  gsl_vector_mul(v,coeff);
  return gsl_blas_dasum(v);
}

/* Returns variance of mix of exponential distributions */
double varExpMix(const gsl_vector* coeff, const gsl_vector* lambda) {
  gsl_vector *v=gsl_vector_alloc(lambda->size);
  gsl_vector *tc=timeConstCube(lambda);
  double mean=meanExpMix(coeff, lambda), EtSqr;
  gsl_vector_memcpy(v, coeff);
  gsl_vector_add(v, coeff); /*coeffs times 2 */
  gsl_vector_mul(v,tc); /* times timeconstants*/
  EtSqr=gsl_blas_dasum(v);
  return EtSqr-(mean*mean);
}

/* Returns standard deviation of mix of exponential distributions */
double stdExpMix(const gsl_vector* coeff, const gsl_vector* lambda) {
  return sqrt(varExpMix(coeff, lambda));
}

void processArgs(int argc, char **argv){
  int i,j;
  int nArgs=4, nArgsMax=4;
  #ifdef OPEN
  const char *usage="modelfile nStates nOpen";
  #elif CLOSED
  const char *usage="modelfile nStates nClosed";
  #endif
  char *modelFn;
  char * endptr;
  FILE *fp;
  
  if( (argc<nArgs) || (argc>nArgsMax)) {
    fprintf(ERR, "Usage:\n%s %s\n", argv[0], usage);
    exit(1);
  }

  modelFn=argv[1];

  nStates=strtol(argv[2], &endptr, 10);
  if(*endptr != '\0') {
    fprintf( ERR, "Error converting %s to integer\n", argv[2]);
    exit(1);
  }

  nDist=strtol(argv[3], &endptr, 10);
  if(*endptr != '\0') {
    fprintf( ERR, "Error converting %s to integer\n", argv[3]);
    exit(1);
  }

  double (*qr)[nStates]=malloc(sizeof (double[nStates][nStates]));
  fp=fopen(modelFn,"r");
  
  io_readMatrix(fp, nStates, nStates, qr);
  /* copy to Qrates */
  Q=gsl_matrix_alloc(nStates, nStates);
  for(i=0; i<nStates; i++) {
    for(j=0; j<nStates; j++) {
      gsl_matrix_set(Q,i,j,qr[i][j]);
    }
  }
  makeConservative(Q);
  fclose(fp);
  free(qr);
}


int main(int argc, char **argv) {
  gsl_vector *lambda, *coeff;

  processArgs(argc, argv);
  
  lambda=gsl_vector_alloc(nDist), coeff=gsl_vector_alloc(nDist);

#ifdef CLOSED
  closedDistribution(Q,lambda, coeff);
#elif OPEN
  openDistribution(Q,lambda, coeff);
#endif

  /* if(showMean) {fprintf(OUT, "%g\n", meanExpMix(coeff, lambda));} */
  /* else if(showVar) {fprintf(OUT, "%g\n", varExpMix(coeff, lambda));} */
  /* else if(showStd) {fprintf(OUT, "%g\n", stdExpMix(coeff, lambda)); return 0;} */
  /* else */
  
  printDist(coeff, lambda);

  return 0;  
}
 
