#ifndef TOOLS
#define TOOLS
/*
 * mctools.h
 *
 *
 * Tools for markov chains
 * 
 * 
 *
 * Ivo Siekmann, 28/04/2009
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
/* #include <gsl/gsl_blas.h> */
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>

/*
 * Fill diagonal so that Q is the matrix of a conservative matrix
 * model.
 */
void makeConservative(gsl_matrix *Q);

/* Compute equilibrium distribution of a markov model */
gsl_matrix* statMarkov(const gsl_matrix *Q);

/* Compute equilibrium distribution of a markov model */
void statMarkovVector(const gsl_matrix *Q, gsl_vector *p);

/* logarithm of normal distribution for x WITHOUT SCALING */
double logNorm1(double x, double mean, double var);

/*Logarithm of the normal distribution WITHOUT SCALING*/
double logNormal(double *q[], const double qMean[], 
		 const double qVar[], int n);

/* Compute the closed times probability distribution. It is assumed
 * that the closed states start at index 0 and are as many as the size
 * of ev.The eigenvalues are saved in ev, the coefficients of the
 * exponentials are saved in coeff.
 */
void closedDistribution(const gsl_matrix* Q, 
		      gsl_vector *ev, gsl_vector *coeff);

/* Compute the open times probability distribution. It is assumed
 * that the closed states start at index 0 and are as many as the size
 * of ev.The eigenvalues are saved in ev, the coefficients of the
 * exponentials are saved in coeff.
 */
void openDistribution(const gsl_matrix* Q, 
		      gsl_vector *ev, gsl_vector *coeff);

/* Computes aggregated distribution from qAA, qAB, qBA and pB. The
 * eigenvalues are saved in ev, the coefficients of the exponentials
 * are saved in coeff. */
void aggregatedStates(const gsl_matrix* qAA, const gsl_matrix* qAB, 
		const gsl_matrix* qBA, const gsl_vector* pB,
		      gsl_vector* ev, gsl_vector *coeff, int Trans);

#endif
