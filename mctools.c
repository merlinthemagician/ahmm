/*
 *  mctools.c
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

#include <stdio.h>
#include <math.h>
/* #include <float.h> */

#ifndef HPC
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_rng.h>
#else
#include "gsl/gsl_complex.h"
#include "gsl/gsl_complex_math.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_permutation.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_rng.h"
#endif

#include "mctools.h"
#include "matrixExp.h"
#include "matrixArith.h"

/* #ifdef TEST */
#include "matrixIO.h"
/* #endif */

#define OUT stdout
#define ERR stderr

#define EPS 1e-12

static gsl_permutation *permut;

static gsl_eigen_nonsymmv_workspace *eigenW;

static gsl_matrix *Acopy;

static gsl_vector_complex *lambda;

static gsl_matrix_complex *eigenT, *eigenTinv;

/* initialise routines where eigenvalues are involved */
void initEigen(int n) {
  eigenW=gsl_eigen_nonsymmv_alloc(n);
  Acopy=gsl_matrix_alloc(n,n);
  lambda=gsl_vector_complex_alloc(n);
  eigenT=gsl_matrix_complex_alloc(n,n);
  eigenTinv=gsl_matrix_complex_alloc(n,n);
  permut=gsl_permutation_alloc(n);
}

/* initialise routines where eigenvalues are involved */
void freeEigen() {
  gsl_eigen_nonsymmv_free(eigenW);
  gsl_matrix_free(Acopy);
  gsl_vector_complex_free(lambda);
  gsl_matrix_complex_free(eigenT);
  gsl_matrix_complex_free(eigenTinv);
}

/* Copies real part of vector cv in v - Error checking*/
void copyVectorReal(gsl_vector *v, const gsl_vector_complex *cv) {
  int i;
  gsl_complex z;

  for(i=0; i<cv->size; i++) {
    z=gsl_vector_complex_get(cv,i);
    if(fabs(GSL_IMAG(z)) > EPS) {
      fprintf(ERR, "warning: mctools: copyVectorReal(): ");
      fprintf(ERR, "copied Vector is not real!");
      exit(1);
    }

    gsl_vector_set(v, i, GSL_REAL(z));
  }
}

void copyMatrixReal(gsl_matrix *M, const gsl_matrix_complex *cM) {
  int i, j;
  gsl_complex z;

  for(i=0; i<cM->size1; i++) {
    for(j=0; j<cM->size2; j++) {
      z=gsl_matrix_complex_get(cM,i, j);
      gsl_matrix_set(M, i, j, GSL_REAL(z));
    }
  }
}

/* Eigensystem computation assuming real eigenvalues and -vectors

   Assumptions: space for n reserved with initEigen().
 */
void realEigensystem(const gsl_matrix *A, gsl_vector *lam, gsl_matrix *V) {
  gsl_matrix_memcpy(Acopy, A);

  gsl_eigen_nonsymmv(Acopy, lambda, eigenT, eigenW);
  copyVectorReal(lam, lambda);
  copyMatrixReal(V,eigenT);
}

/* inverts A using the LU transform */
void invertMatrix(const gsl_matrix *A, gsl_matrix *Ainv) {
  int sig;
  gsl_matrix_memcpy(Acopy, A);

  gsl_linalg_LU_decomp(Acopy, permut, &sig);

  gsl_linalg_LU_invert(Acopy, permut, Ainv);
}

/* computes the i-th summand of the spectral expansion if the trafo
   matrices are T and Tinv */
void spectralExpansion(const gsl_matrix *T, const gsl_matrix *Tinv, int i,
		       gsl_matrix *A) {
  gsl_vector_const_view cv = gsl_matrix_const_column(T, i);

  /* Compute T Proj[i] Tinv */
  gsl_matrix_set_zero(Acopy);
  gsl_matrix_set_col(Acopy, i, &cv.vector);
  
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Acopy, Tinv, 0.0, A);
}


/* Computes aggregated distribution from qAA, qAB, qBA and pB. The
 * eigenvalues are saved in ev, the coefficients of the exponentials
 * are saved in coeff. */
void aggregatedStates(const gsl_matrix* qAA, const gsl_matrix* qAB, 
		      const gsl_matrix* qBA, const gsl_vector* pB,
		      gsl_vector* ev, gsl_vector *coeff, int Trans) {

  int i, N=ev->size;
  gsl_matrix *T, *Tinv, *Ai, *gABi;
  gsl_vector *phi;
  gsl_vector *lambda, *pVecI;
  double c;

/*   fprintf(OUT, "aggregated(): %i aggregated states.\n", N); */
  phi = gsl_vector_alloc(N);
  lambda=gsl_vector_alloc(N);
  pVecI=gsl_vector_alloc(N);
  
  T=gsl_matrix_alloc(N,N);  
  Tinv=gsl_matrix_alloc(N,N);
  Ai=gsl_matrix_alloc(N,N);
  /* Size of B is qAB->size2 ?!? */
  gABi=gsl_matrix_alloc(N, (qAA->size1)+(pB->size)-N);

/*   fprintf(OUT, "Initialising done.\n"); */
  
  initEigen(N); 

/*   fprintf(OUT, "qBA (%i x %i), pB (%i), phi (%i)\n", */
/* 	  qBA->size1,qBA->size2,pB->size, phi->size); */
  gsl_blas_dgemv (CblasTrans,1.0,qBA,pB,0,phi);

  toProbVector(phi);

  realEigensystem(qAA, lambda, T);
/*   fprintf(OUT, "Eigensystem computed:\n"); */
/*   for(i=0; i<lambda->size; i++) fprintf(OUT, "\t%f",gsl_vector_get(lambda, i)); */
/*   fprintf(OUT, "\n"); */

  invertMatrix(T, Tinv);
/*   fprintf(OUT, "Matrix inverted\n"); */

  for(i=0; i<lambda->size; i++) {
    spectralExpansion(T, Tinv, i, Ai);

/*     fprintf(OUT, "Ai (%i x %i), qAB (%i x %i), gABi (%i x %i)\n", */
/* 	    Ai->size1, Ai->size2, qAB->size1, qAB->size2, */
/* 	    gABi->size1, gABi->size2 */
/* 	    ); */

    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Ai, qAB, 0.0, gABi);

    /*SPALTEN VON gABi SUMMIEREN! --> (Colqhoun & Hawkes, 1983; eq. 1.24 )*/
/*     fprintf(OUT, "Sum columns of gABi\n"); */
/*     printMat(OUT, gABi); */

    sumColumns(gABi, pVecI);

    /* GEWICHTET MIT MOEGLICHEN STARTZUSTAENDEN SUMMIEREN */
    gsl_blas_ddot(pVecI, phi, &c);
    gsl_vector_set(coeff, i, c);
    
  }

  gsl_vector_memcpy(ev, lambda);

  freeEigen();
}

/* Compute the closed times probability distribution. It is assumed
 * that the closed states start at index 0 and are as many as the size
 * of ev.The eigenvalues are saved in ev, the coefficients of the
 * exponentials are saved in coeff.
 */
void closedDistribution(const gsl_matrix* Q, 
			gsl_vector *ev, gsl_vector *coeff) {
  int N = ev->size;

  gsl_matrix_const_view qAA = 
    gsl_matrix_const_submatrix(Q, 0, 0, N, N);
  gsl_matrix_const_view qAB = 
    gsl_matrix_const_submatrix(Q, 0, N, N, Q->size2-N);
  gsl_matrix_const_view qBA =
    gsl_matrix_const_submatrix(Q, N, 0, Q->size1-N, N);

  gsl_vector *statD=gsl_vector_alloc(Q->size1);
  gsl_vector_const_view pB=gsl_vector_const_subvector(statD,N,Q->size1-N);

  statMarkovVector(Q, statD);

/*   fprintf(OUT, "Compute closed distribution.\n"); */
  aggregatedStates(&qAA.matrix, &qAB.matrix, &qBA.matrix, &pB.vector, 
		   ev, coeff, 1);
}

/* Compute the open times probability distribution. It is assumed
 * that the closed states start at index 0 and are as many as the size
 * of ev.The eigenvalues are saved in ev, the coefficients of the
 * exponentials are saved in coeff.
 */
void openDistribution(const gsl_matrix* Q, 
		      gsl_vector *ev, gsl_vector *coeff) {
  int N = ev->size;

  /* qBB */
  gsl_matrix_const_view qAA = 
    gsl_matrix_const_submatrix(Q, Q->size1-N, Q->size1-N, N, N);
  /*qBA*/
  gsl_matrix_const_view qAB = gsl_matrix_const_submatrix(Q, Q->size1-N, 0, N, Q->size2-N);
  /*qAB*/
  gsl_matrix_const_view qBA =gsl_matrix_const_submatrix(Q, 0, Q->size2-N, Q->size1-N, N);

  gsl_vector *statD=gsl_vector_alloc(Q->size1);
  gsl_vector_const_view pB=gsl_vector_const_subvector(statD,0,Q->size1-N);

  statMarkovVector(Q, statD);
/*   fprintf(OUT, "Compute open distribution.\n"); */
  aggregatedStates(&qAA.matrix, &qAB.matrix, &qBA.matrix, &pB.vector, 
		   ev, coeff, 1);
}


/* null space of a square matrix --- not implemeted yet! */
void nullSpaceLU(gsl_matrix* M) {
  gsl_permutation *p=gsl_permutation_alloc(M->size1);
  int sig;

  gsl_linalg_LU_decomp(M, p, &sig);

/*   printMat(OUT,M); */

}

/* returns the null space of a square matrix as a (Corang x N) matrix*/
gsl_matrix * nullSpaceSVD(gsl_matrix* M) {
  double eps=1e-12;
  int d=M->size2;
  int corg;
  gsl_vector *s=gsl_vector_alloc(d), *w=gsl_vector_alloc(d);
  gsl_matrix *V=gsl_matrix_alloc(d, d);
  gsl_matrix_view singView;
  gsl_matrix *ns;

/*   fprintf(OUT, "Dimension of M: %i\n", d); */
/*   printMat(OUT, M); */

  gsl_linalg_SV_decomp(M, V,s, w);
/*   fprintf(OUT, "Decomposition done.\n"); */

  for(corg=0; (corg<d) && (gsl_vector_get(s,d-corg-1)<eps); corg++);
/*   vprint(OUT, s); */
/*   fprintf(OUT,"\nCo-rank: %i\n",corg); */

  ns=gsl_matrix_alloc(corg,d);
  singView=gsl_matrix_submatrix(V,0,d-corg, M->size1, corg);
  gsl_matrix_transpose_memcpy(ns,&singView.matrix);  

/*   printMat(OUT,V); */
/*   fprintf(OUT, "\n"); */
/*   printMat(OUT,ns); */

  gsl_vector_free(s), gsl_vector_free(w);
  gsl_matrix_free(V);

  return ns;

}

/*normalise columns assuming that Q is conservative*/
void normCols(gsl_matrix *Q) {
  int i,j;
  double sum;
  for(j=0; j<Q->size1; j++) {
    sum=-gsl_matrix_get(Q,j,j);
    for(i=0;i<Q->size2;i++) {
      gsl_matrix_set(Q,i,j,gsl_matrix_get(Q,i,j)/sum);
    }
  }
}

/*normalise columns assuming that Q is conservative*/
void normRows(gsl_matrix *Q) {
  int i,j;
  double sum;
  for(i=0; i<Q->size1; i++) {
    sum=-gsl_matrix_get(Q,i,i);
    for(j=0;j<Q->size2;j++) {
      gsl_matrix_set(Q,i,j,gsl_matrix_get(Q,i,j)/sum);
    }
  }
}

/* Compute equilibrium distribution of a markov model */
void statMarkovVector(const gsl_matrix *Q, gsl_vector *p) {
  gsl_matrix* trQ= gsl_matrix_alloc(Q->size1, Q->size2);
  gsl_matrix *res;

  gsl_matrix_transpose_memcpy (trQ,Q);
  res=nullSpaceSVD(trQ);
  
  if(res->size1 > 1) {
    fprintf(ERR,"No unique equilibrium distribution.\n");
    exit(1);
  }
  
  /* fprintf(OUT, "Computing stationary distribution of:\n"); */
  /* printMat(OUT, Q); */
  gsl_matrix_get_row(p,res, 0);
  toProbVector(p);
  gsl_matrix_free(res), gsl_matrix_free(trQ);

}

/* Compute equilibrium distribution of a markov model */
gsl_matrix* statMarkov(const gsl_matrix *Q) {
  int i;
  gsl_matrix* trQ= gsl_matrix_alloc(Q->size1, Q->size2);
  gsl_matrix*res;
  gsl_vector_view rowView;

  gsl_matrix_transpose_memcpy (trQ,Q);

  res=nullSpaceSVD(trQ);

  for(i=0; i<res->size1; i++) {
    rowView=gsl_matrix_row(res, i);
    toProbVector(&rowView.vector);
  }

  return res;
}

/*Fill diagonal so that Q is a matrix of a conservative matrix model.
 */
void makeConservative(gsl_matrix *Q) {
  int i, j;
  double sumOfRates = 0;
  for(i=0; i<Q->size1; i++) {
    sumOfRates = 0;
    for(j=0; j<i; j++) {
      sumOfRates += gsl_matrix_get(Q,i,j);
    }
    for(j=i+1; j<Q->size2; j++) {
      sumOfRates += gsl_matrix_get(Q,i,j);
    }
    gsl_matrix_set(Q,i,i, -sumOfRates);
  }
}

/* logarithm of normal distribution for x WITHOUT SCALING */
double logNorm1(double x, double mean, double var) {
  return - (x-mean)*(x-mean)/(2*var);
}

/*Logarithm of the normal distribution WITHOUT SCALING*/
double logNormal(double *q[], 
		 const double qMean[], const double qVar[], int n) {
  double sum=0;
  int i;

  for(i=0; i<n; i++) {
    sum += logNorm1(q[i][0], qMean[i], qVar[i]);
  }
  return sum;
}


