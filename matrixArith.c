/*
 *  matrixArith.c
 *
 *
 * Matrix operation
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

#include "matrixIO.h"

#define ERR stderr

#define EPS 1e-12

#define FLIP(x) ((x)==0)?1:0

/* Produkt der Matrizen M1 und M2:
 * R darf NICHT mit einer der beiden Matrizen M1 oder M2
 * uebereinstimmen.
 */
void multMat(const gsl_matrix *M1, const gsl_matrix *M2, gsl_matrix *R) {
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
		 1.0, M1, M2, 0.0, R);
}

/* Quadrat der Matrix M
 * R darf NICHT mit M uebereinstimmen.
 */
void sqrMat(const gsl_matrix *M, gsl_matrix *R) {
  multMat(M, M, R);
}


/* Frobenius norm of matrix M */
double frobeniusNorm(const gsl_matrix *M) {
  double sqrsum=0, mIJ;
  int i, j;

  for(i=0; i<M->size1; i++) {
    for(j=0; j<M->size2; j++) {
      mIJ=gsl_matrix_get(M, i, j);
      sqrsum += mIJ*mIJ;
    }
  }

  return sqrt(sqrsum);
}


/* Efficient matrix power function, implemented by repeated squaring */ 
void matrixPower(const gsl_matrix * A, int n, gsl_matrix *An) {
  /*refering to which of the two B?*/
  int lr=0;

  /*An or res ?*/
  int toAn=0;

  gsl_matrix *B[2];
  gsl_matrix *res=gsl_matrix_alloc(A->size1, A->size2);

  gsl_matrix_set_identity (An);

  if( n== 0) {
    /*An --> IdentityMatrix */
    return;
  }

/*   gsl_matrix_memcpy(An, A); */

  B[0] = gsl_matrix_alloc(A->size1, A->size2);
  B[1] = gsl_matrix_alloc(A->size1, A->size2);
  
  gsl_matrix_memcpy(B[lr], A);

  while(n > 1) {
    if(n%2 == 0) {
      n /=2;
      sqrMat(B[lr], B[FLIP(lr)]);
      lr=FLIP(lr);
/*       fprintf(OUT,"\nSqrM:\n"); */
/*       printMat(OUT, B[lr]); */
    }
    else {
      n--;
      
      if (toAn){
	multMat(res, B[lr], An);
	toAn=FLIP(toAn);
/* 	fprintf(OUT,"\nres:\n"); */
/* 	printMat(OUT, res); */
      }
      else {
	multMat(An, B[lr], res);
	toAn=FLIP(toAn);
/* 	fprintf(OUT,"\nAn:\n"); */
/* 	printMat(OUT, An); */
      }
    }
  }
  /*if current result in res: copy to An*/
  if(toAn){
    multMat(res, B[lr], An);
  }
  else {
    multMat(An, B[lr], res);
    gsl_matrix_memcpy(An,res);
  }
    /* gsl_matrix_memcpy(An, B[lr]); */
  gsl_matrix_free(B[0]),  gsl_matrix_free(B[1]), gsl_matrix_free(res);
}

/* sum columns of a matrix -> vector*/
void sumColumns(const gsl_matrix *M, gsl_vector *v) {
  int j;
  gsl_matrix_get_col(v, M, 0);

  for(j=1; j<M->size2; j++) {
    gsl_vector_const_view cv=gsl_matrix_const_column(M, j);
    gsl_vector_add(v, &cv.vector);
  }
}

/* scales v so that the sum of its components is 1 - returns scaling factor */
double toProbVector(gsl_vector *v) {
  int i;
  double sum=0;

  if( (! gsl_vector_isnonneg(v)) && (! gsl_vector_isneg(v))) {
    fprintf(ERR, "toProbVector(): Scaling not possible:\n");
    vprint(ERR,v);
    fprintf(ERR, "Vector is neither nonnegative nor negative!\n");
    exit(1);
  }

  for(i=0; i<v->size; i++)
    sum += gsl_vector_get(v,i);
  if(fabs(sum) <= EPS) {
    fprintf(ERR, "toProbVector(): cannot divide by %f\n", sum);
    exit(1);
  }
  gsl_vector_scale(v,1.0/sum);

  return sum;
}


/*Set n rows from bottom in matrix M to zero, save in PM*/
void zeroRowsBottom(const gsl_matrix *M, int n, gsl_matrix *PM) {
  int i;
  gsl_vector_view row;

  gsl_matrix_memcpy(PM, M);
  for(i=0; i<n; i++) {
    row = gsl_matrix_row(PM, PM->size1-1-i);
    gsl_vector_set_zero(&row.vector);
  }
}

/*Set n rows from top in matrix M to zero, save in PM*/
void zeroRowsTop(const gsl_matrix *M, int n, gsl_matrix *PM) {
  int i;
  gsl_vector_view row;

  gsl_matrix_memcpy(PM, M);
  for(i=0; i<n; i++) {
    row = gsl_matrix_row(PM, i);
    gsl_vector_set_zero(&row.vector);
  }
}

/* relative error with sign */
void relativeErrorSigned(const gsl_matrix *M, 
			 const gsl_matrix *Mest, gsl_matrix *res) {
  gsl_matrix_memcpy(res, Mest);
  gsl_matrix_sub(res, M);
  gsl_matrix_div_elements(res,M);
}

/* relative error with sign in percent */
void relativeErrorSignedPercent(const gsl_matrix *M, 
				const gsl_matrix *Mest, gsl_matrix *res) {

  relativeErrorSigned(M, Mest, res);
  gsl_matrix_scale(res,100);
}

