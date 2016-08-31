#ifndef HMMUTILS
#define HMMUTILS
/*
 * mp_HMMutils.c
 *
 * Helper routines for HMM
 * 
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

#include <strings.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>

#include "mp_parameter.h"

#include "likelihood.h"

#define NDATA 10000000

/* Exponential prior */
double expPrior (const parameters *p, int nP,
		 const intparameters *ip, int nIP);

/* Exponential prior */
double expDoublePrior (const parameters *p, int nP, likelihood *L);

/* Reads last line of a text file */
char *lastLine(FILE *fp, char *buf, size_t max_len);

/* Convert s containing ints to an int vector*/
size_t str2ints(char* s, int *nums);

/* Convert s containing double to a double vector*/
size_t str2doubleV(char* s, double *nums);

/*Reads a list of doubles*/
int readDoubles(FILE *fp, double *v, int buf);

/*Reads a dynamically allocated double matrix. Returns number of
  rows. */
int readMatrix(FILE *fp, double **A, int M, int N);

/*Generate O/C trace from thresholded open probabilities: Po is in column 0, length of segment is in column 2 */
void modes2Events(int *events, const double **A, double Pthreshold,
		  int nA);

/* Estimates the open sequences and the length of OCO sequences by
   thresholding and saves to NC and NO respectively - returns the
   length of the NC "histogram" */
int guessOClength(const double *data, int *NC, int *NO, int *n, double thresh, double *pC);

/* Reads data from a file datafp with nRows rows. Threshold must be
   provided in thresh.  Returns the NC "histogram" and saves the
   number of open events to NO. Maximum number of closed events saved
   to maxClosed, ACTUAL number of rows read passed to nRows, total
   size of data set saved in nData, estimated closed probability saved
   in pC.*/
int* columnsData2NOC(const char *datafn, double thresh, int *NO, int *maxClosed,
		     int* nRows, int *nData, double *pC);

/* Reads data from a file datafp with nRows rows. Threshold must be
   provided in thresh.  Returns events, an integer array of zeroes of
   ones. The ACTUAL number of rows read passed to nRows, total size of
   data set saved in nData.*/
int* mp_columnsData2Events(const char *datafn, double thresh,
			   int* nRows, int *nData);

int *mp_matrix2events(FILE* datafp, int *nRows, int nCol, int *nEvents, double Pthresh);

/* Converts A to NO and NC histogram: Po is in column 0 of A, 
   length of segment is in column 2. */
int modes2OClength(const double **A, int *NC, int *NO, int mA, double Pthresh, double *pC);

int* mp_matrix2OClength(FILE* datafp, int *nRows, int nCol,
			int *NO, int *maxClosed,
			double Pthresh, double *pC);

/* Read data from file fp, threshold values for obtaining a sequence
   of open and closed events. Calculate sequences of closed events NC
   and calculate number of open states NO. Calculate closed
   probability pC.Return maximum length of closed sequence. */
int processDataOneOpen(FILE *fp, int *NC, int *NO, int *nData, double thresh, double *pC);


/* count positive entries in Q */
int countPositive(const gsl_matrix *Q);

/* converts Markov model to parameters: qIJ in first half, qJI in second half */
/* Q cannot be const because org[k] pointers are not const */
void matrix2ParameterSubdiag(gsl_matrix *Q, gsl_matrix *Qp, 
			     parameters *p, int nPar,int *qInd);

/* converts Markov model to parameters: qIJ at even, qJI in odd
   indices. Q cannot be const because org[k] pointers are not const */ 
void matrix2ParameterPairs(gsl_matrix *Q, gsl_matrix *Qp, 
			   parameters *p, int nPar,int *qInd);

/* converts Markov model to parameters: qIJ is saved at even, qJI at odd
   indices starting from k0 */
void matrix2ParameterPairsK0(gsl_matrix *Q, gsl_matrix *Qp, 
			     parameters *p, 
			     int k0, int nPar,int *qInd);

/*output rates and additionally calculate stationary distribution*/
void QandPhiout(FILE *ratesfp, int iteration, 
		const parameters *rates, int nP, 
		FILE *statfp, const gsl_matrix *Q);

/* Gibbs sampling strategy for sampling rates and rate constants */
/* alpha is the set of hyperparameters for the rate constants */
/* ONE COULD GIVE a vector which is a proposal for the stationary
   distribution as an argument */
void samplePhiandQ(const gsl_rng *r, 
		   parameters *p, int nPar,
		   const double *alpha,
		   gsl_matrix *Q,
		   const int *qIndices);

/* Sample steps parameter according to constraints starting from n0 */
void sampleN(const gsl_rng *r, 
	     parameters *p, int n0, int steps);

/* Enforces detailed balance for a given matrix with cycle */
void forceDetailedBalance(gsl_matrix *Qp,
			  const int *cycle, int nCycle);

/* Checks if a model Q satisfies the detailed balance condition */
int checkDetailedBalance(const gsl_matrix *Q);


/* Read model from file modelfile in double vector model */
double * readModel(const char *modelfile, int nStates);

/* Get Arguments */
void mp_getArgsOneOpen(char **argv, int argc,char **datafn,
		       char **modelfile,
		       int *nStates,
		       int *nIter,
		       double *delta,
		       int *seed,
		       char **restartFn,
		       double *thresh,
		       char **ratesFn,
		       char **statFn,
		       char **likeFn);

/* Get Arguments for N open states. */
void mp_getArgsNOpen(char **argv, int argc,char **datafn,
		     char **modelfile,
		     int *nStates,
		     int *nOpen,
		     int *nIter,
		     double *delta,
		     int *seed,
		     char **restartFn,
		     double *samplingInt,		     
		     double *thresh,
		     char **ratesFn,
		     char **statFn,
		     char **likeFn);
#endif
