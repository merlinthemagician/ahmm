#ifndef ARITH
#define ARITH
/*
 *  matrixArith.c
 *
 *
 * Matrix operations
 * 
 * 
 *
 * Ivo Siekmann, 21/07/2009
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */

#include <stdio.h>
#include <math.h>

#ifndef HPC
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#else
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_blas.h"
#endif

/* Produkt der Matrizen M1 und M2:
 * R darf NICHT mit einer der beiden Matrizen M1 oder M2
 * uebereinstimmen.
 */
void multMat(const gsl_matrix *M1, const gsl_matrix *M2, gsl_matrix *R);

/* Quadrat der Matrix M
 * R darf NICHT mit M uebereinstimmen.
 */
void sqrMat(const gsl_matrix *M, gsl_matrix *R);


/* Frobenius norm of matrix M */
double frobeniusNorm(const gsl_matrix *M);

/* Efficient matrix power function, implemented by repeated squaring */ 
void matrixPower(const gsl_matrix * A, int n, gsl_matrix *An);

/* sum columns of a matrix -> vector*/
void sumColumns(const gsl_matrix *M, gsl_vector *v);

/* scales v so that the sum of its components is 1 */
double toProbVector(gsl_vector *v);

/*Set n rows from bottom in matrix M to zero, save in PM */
void zeroRowsBottom(const gsl_matrix *M, int n, gsl_matrix *PM);

/*Set n rows from top in matrix M to zero, save in PM*/
void zeroRowsTop(const gsl_matrix *M, int n, gsl_matrix *PM);

/* relative error with sign */
void relativeErrorSigned(const gsl_matrix *M, 
			 const gsl_matrix *Mest, gsl_matrix *res);

/* relative error with sign in percent */
void relativeErrorSignedPercent(const gsl_matrix *M, 
				const gsl_matrix *Mest, gsl_matrix *res);
#endif
