#ifndef HMMlikelihoods
#define HMMlikelihoods
/*
 * mp_parameter.c
 *
 * MCMC parameter data type
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
#include <gsl/gsl_matrix.h>

/* Computes the likelihood difference of Q and Qp for a given
 * sequence of consecutive closed and open states. */
double logprNOCdiff1(const int NO, const int *NC, 
		     gsl_matrix *Q, 
		     gsl_matrix *Qp, double tau, int nC);
/* Computes the likelihoods of Q and Qp for a given
 * sequence of consecutive closed and open states. */
double logprNOC1(int NO, int *NC, 
		     gsl_matrix *Q, 
		 gsl_matrix *Qp, double tau, int nC, likelihood *L);
#endif
