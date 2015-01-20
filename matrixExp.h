#ifndef EXP
#define EXP
/*
 *  matrixExp.h
 *
 *  Pade approximation with scaling and squaring
 *
 * Ivo Siekmann, 21/07/2009
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */

#include <gsl/gsl_matrix.h>

/* "diagonal" Pade approximation of matrix exponential and scaling and
   squaring using parameters according to (Moler & van Loan, 2003)*/
void padeExp(const gsl_matrix *A, double tau, gsl_matrix *expA);

#endif
