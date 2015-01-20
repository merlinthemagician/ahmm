/*
 *  matrixExp.c
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

#include <stdio.h>
#include <math.h>

#ifndef HPC
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#else
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_permutation.h"
#include "gsl/gsl_linalg.h"
#endif

#include "matrixIO.h"
#include "matrixArith.h"

#define OUT stdout
#define ERR stderr

/* Computes binomial(n, k) recursively. */
int binomial(int n, int k) {
  if( (k==0) || (k==n) || (n==1) || (n==0) ) return 1;
  if( (k==1) || (k==n-1)) return n;

  return binomial(n-1, k) + binomial(n-1, k-1);
}

/* Computes jth pade coefficient of order q, consists of a binomial
   coefficient divided by a factorial */
double padeCoeff(int q, int k) {
  double denom=1, mul=2*q;
  int i;

  if(k==0) return 1;
  if(k==1) return 0.5;

  for(i=1, mul=2*q; i<=k; i++, mul--) {
    denom*=mul;
  }
  return binomial(q,k)/denom;
}


/* "diagonal" Pade approximation of matrix exponential, order q */
void diagonalPade(int q, const gsl_matrix *A, gsl_matrix *expA) {
  gsl_matrix *D = gsl_matrix_alloc(A->size1, A->size2);
  gsl_matrix *M = gsl_matrix_alloc(A->size1, A->size2);
  gsl_matrix *Offs=gsl_matrix_alloc(A->size1, A->size2);
  gsl_matrix *An=gsl_matrix_alloc(A->size1, A->size2);
  int i;
  gsl_permutation  *p = gsl_permutation_alloc(A->size1);

  gsl_matrix_set_identity(D), gsl_matrix_set_identity(M);
  gsl_matrix_set_identity(An);
/*   gsl_matrix_memcpy(Offs, A); */
  
  for(i=1; i<=q; i++) {
    multMat(An, A, Offs);
    gsl_matrix_memcpy(An, Offs);
    gsl_matrix_scale(Offs, padeCoeff(q, i));

    gsl_matrix_add(M, Offs);
/*     fprintf(OUT, "M+Offs:\n"); */
/*     printMat(OUT, M); */
/*     fprintf(OUT, "\n"); */

    if( (i%2) == 1 ) gsl_matrix_sub(D, Offs); /* gsl_matrix_scale(Offs, -1); */
    else gsl_matrix_add(D, Offs);

/*     fprintf(OUT, "D+(-1)^j (Offs):\n"); */
/*     printMat(OUT, D); */
/*     fprintf(OUT, "\n"); */
  }
/*   fprintf(OUT, "M:\n"); */
/*   printMat(OUT, M); */
/*   fprintf(OUT, "\n"); */

/*   fprintf(OUT, "D:\n"); */
/*   printMat(OUT, D); */
/*   fprintf(OUT, "\n"); */

  gsl_linalg_LU_decomp (D, p, &i);
/*   fprintf(OUT, "RowEchelon D:\n"); */
/*   printMat(OUT, D); */
/*   fprintf(OUT, "\n"); */

  gsl_linalg_LU_invert(D, p, Offs);
/*   fprintf(OUT, "RowEchelon D:\n"); */
/*   printMat(OUT, Offs); */
/*   fprintf(OUT, "\n"); */

  multMat(Offs, M, expA);

  gsl_matrix_free(D), gsl_matrix_free(M);
  gsl_matrix_free(An), gsl_matrix_free(Offs);
  gsl_permutation_free(p);
}

/* "diagonal" Pade approximation of matrix exponential, order q,
   scaling and squaring k times */
void dPadeSclSqr(int q, int k, const gsl_matrix *A, gsl_matrix *expA) {
  gsl_matrix *Acpy=gsl_matrix_alloc(A->size1, A->size2);
  gsl_matrix *expAcpy = gsl_matrix_alloc(A->size1, A->size2);
  double scl=pow(0.5, k);
  int i;

  if(k <= 0) {
    diagonalPade(q,A,expA);
    return;
  }

  gsl_matrix_memcpy(Acpy, A);
  gsl_matrix_scale(Acpy, scl);

  diagonalPade(q, Acpy, expAcpy);

  for(i=1;i<=k;i++) {
    if( i%2==0 ) {
      sqrMat(Acpy,expAcpy);
/*       printMat(OUT, expAcpy); */
/*       fprintf(OUT, "\n"); */
    }
    else {
      sqrMat(expAcpy, Acpy);
/*       printMat(OUT, Acpy); */
/*       fprintf(OUT, "\n"); */
    }
  }

  if( k%2 == 0 )
    gsl_matrix_memcpy(expA, expAcpy);
  else 
    gsl_matrix_memcpy(expA, Acpy);

  gsl_matrix_free(expAcpy), gsl_matrix_free(Acpy);
}

static int qm[6][2] = {{3,0},{4,0},{6,1},{6,5},{6,8},{6,11}};

/* "diagonal" Pade approximation of matrix exponential and scaling and
   squaring using parameters according to (Moler & van Loan, 2003)*/
void padeExp(const gsl_matrix *A, double tau, gsl_matrix *expA) {
  gsl_matrix *Acpy = gsl_matrix_alloc(A->size1, A->size2);
  double norm;

  gsl_matrix_memcpy(Acpy, A);
  gsl_matrix_scale(Acpy, tau);
  norm=frobeniusNorm(Acpy);
/*   fprintf(OUT, "Matrix norm (%g)", norm); */
  if(norm<=1e-2){
/*     fprintf(OUT, " below 1e-2: (%i,%i)\n", qm[0][0], qm[0][1]); */
    dPadeSclSqr(qm[0][0], qm[0][1], Acpy, expA);
  }
  else if( (norm > 1e-2) && (norm <=1e-1)) {
/*     fprintf(OUT, " between 1e-2 and 0.1: (%i,%i)\n", qm[1][0], qm[1][1]); */
    dPadeSclSqr(qm[1][0], qm[1][1], Acpy, expA);
  }
  else if( (norm > 1e-1) && (norm <=1)) {
/*     fprintf(OUT, " between 0.1 and 1: (%i,%i)\n", qm[2][0], qm[2][1]); */
    dPadeSclSqr(qm[2][0], qm[2][1], Acpy, expA);
  }
  else if( (norm > 1) && (norm <=1e1)) {
/*     fprintf(OUT, " between 1 and 10: (%i,%i)\n", qm[3][0], qm[3][1]); */
    dPadeSclSqr(qm[3][0], qm[3][1], Acpy, expA);
  }
  else if( (norm > 1e1) && (norm <=1e2)) {
/*     fprintf(OUT, " between 10 and 100: (%i,%i)\n", qm[4][0], qm[4][1]); */
    dPadeSclSqr(qm[4][0], qm[4][1], Acpy, expA);
  }
  else if( (norm > 1e2) && (norm <=1e3)) {
/*     fprintf(OUT, " between 100 and 1000: (%i,%i)\n", qm[5][0], qm[5][1]); */
    dPadeSclSqr(qm[5][0], qm[5][1], Acpy, expA);
  }
  else {
    fprintf(ERR, "padeExp(): Pade approximation not supported for ");
    fprintf(ERR, "matrix norms exceeding 1000.\n");
    printMat(OUT,A);
    exit(1);
  }
  gsl_matrix_free(Acpy);
}

#ifdef MAIN
#include "mctools.h"
int main(int argc, char **argv) {
  
  if(argc >= 3) {
    char *modelFN = argv[1];
    FILE *fp=fopen(modelFN, "r");
    int n = atoi(argv[2]);
    double *Q = malloc(n*n*sizeof(double));
    gsl_matrix_view qView = gsl_matrix_view_array(Q, n, n);
    gsl_matrix *expQ = gsl_matrix_alloc(n , n);
    double tau = 1;

    if(argc>=4) {
      tau=atof(argv[3]);
    }

    doubleVec(fp, Q);
    fclose(fp);
    makeConservative(&qView.matrix);

    padeExp(&qView.matrix, tau, expQ);

    fprintf(OUT, "matrix:\n");
    printMat(OUT, &qView.matrix);

    fprintf(OUT, "\nmatrix exponential:\n");
    /* printMat(OUT, expQ); */
    printMatf(OUT, expQ, "\n", "\t", "%g");


    if(argc >= 5) {
      char *modelFN2 = argv[4];
      double *Q2 = malloc(n*n*sizeof(double));
/*       gsl_matrix_view q2View = gsl_matrix_view_array(Q2, n, n); */
/*       gsl_matrix *expQ2 = gsl_matrix_alloc(n , n); */
      gsl_matrix_view expQ2View = gsl_matrix_view_array(Q2, n, n);
      gsl_matrix *err = gsl_matrix_alloc(n , n);

      fp=fopen(modelFN2, "r");
      doubleVec(fp, Q2);

/*       padeExp(&q2View.matrix, tau, expQ2); */

/*       fprintf(OUT, "matrix:\n"); */
/*       printMat(OUT, &q2View.matrix); */

      fprintf(OUT, "\nmatrix exponential:\n");
      printMat(OUT, &expQ2View.matrix);

      gsl_matrix_memcpy(err, &expQ2View.matrix);
      gsl_matrix_sub(err, expQ);
/*       gsl_matrix_scale(err, -1); */
      fprintf(OUT,"\nabsolute difference Q2-Q1:\n");
      printMat(OUT, err);

      fprintf(OUT,"\nRelative error: (Q2-Q1)/Q1\n");
      gsl_matrix_div_elements(err, expQ);
      printMat(OUT, err);

      fprintf(OUT,"\nRelative error: (Q2-Q1)/Q1 (percent)\n");
      gsl_matrix_scale(err, 100);
/*       printMat(OUT, err); */
      printMatf(OUT, err, "\n","\t", "%+.2f");
    }
  }
  else {
    fprintf(ERR, "Usage: matrixExp modelFile N [tau] [model]\n");
    exit(1);
  }

  return 0;
}
#endif


#ifdef TEST
#include <gsl/gsl_rng.h>

int main(int argc, char **argv) {

    int n, k;

  double Q[] = {-0.058,0.058,0,0,
		0.3,-6.9,1.7,4.9,
		0,0.6,-0.6,0,
		0,0.8,0,-0.8};

  gsl_matrix_view qView = gsl_matrix_view_array(Q,4,4);
  gsl_matrix *expQ=gsl_matrix_alloc(4,4);

  for(n=0; n<=6; n++) {
    for(k=0; k<=n; k++) {
      printf("%i ", binomial(n,k));
    }
    printf("\n");
  }
  printf("\n");

  for(n=0; n<=6; n++) {
    for(k=0; k<=n; k++) {
      printf("%g ", padeCoeff(n,k));
    }
    printf("\n");
  }
  printf("\n");



  fprintf(OUT, "Q:\n");
  printMat(OUT, &qView.matrix);
  fprintf(OUT, "\n");

  for(k=2; k<6; k++) {
    fprintf(OUT, "Pade-Approximation %i:\n", k);
    diagonalPade(k,&qView.matrix, expQ);
    printMat(OUT, expQ);
    fprintf(OUT, "\n");
  }

  /* "diagonal" Pade approximation of matrix exponential, order q */
  fprintf(OUT, "Frobenius norm of Q = %f\n", 
	  frobeniusNorm(&qView.matrix));

/*   gsl_matrix_scale(&qView.matrix, 1.0/32); */
/*   fprintf(OUT, "Frobenius norm of scaled Q = %f\n",  */
/* 	  frobeniusNorm(&qView.matrix)); */
 

/* "diagonal" Pade approximation of matrix exponential, order q,
   scaling and squaring k times */
/*   dPadeSclSqr(5, 5, &qView.matrix, expQ); */

/* "diagonal" Pade approximation of matrix exponential and scaling and
   squaring using parameters according to (Moler & van Loan, 2003)*/
  padeExp(&qView.matrix, 1, expQ);

/*   diagonalPade(5,&qView.matrix, expQ); */

/*   sqrMat(expQ, &qView.matrix); */
/*   sqrMat(&qView.matrix, expQ); */
/*   sqrMat(expQ, &qView.matrix); */
/*   sqrMat(&qView.matrix, expQ); */
/*   sqrMat(expQ, &qView.matrix); */
  fprintf(OUT, "Scaling and squaring:\n");
  printMat(OUT, expQ);
  fprintf(OUT, "\n");

  fprintf(OUT, "Random matrix:\n");
  gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
  int i,j;
  for(i=0; i<qView.matrix.size1; i++) {
    for(j=0; j<qView.matrix.size2; j++) {
      double Qij= 20*(gsl_rng_uniform(r)-0.5);
      gsl_matrix_set(&qView.matrix, i, j, Qij);
    }
  }
  printMat(OUT, &qView.matrix);
  fprintf(OUT, "\n");
  fprintf(OUT, "Scaling and squaring:\n");
  padeExp(&qView.matrix, 1, expQ);
  printMat(OUT, expQ);


  return 0;  
}
#endif

