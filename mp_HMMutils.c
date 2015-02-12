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
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef HPC
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#else
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#endif

#include "mp_parameter.h"
#include "mp_proposal.h"
#include "mp_mcmc.h"
#include "likelihood.h"

#include "mctools.h"
#include "matrixIO.h"
#include "nw_data.h"

#define OUT stdout
#define ERR stderr
#define EPS 1e-12

#define NDATA 10000000

static const double expPar=30;

/* Exponential prior */
double expPrior (const parameters *p, int nP,
		 const intparameters *ip, int nIP) {
  double sum=0;
  int i;

  for(i=0; i<nP; i++) {
    double old=getParameter(p,i), new=getProposal(p,i);
    sum+=1/expPar*(-new+old);
  }
  return sum;
}

/* Exponential prior */
double expDoublePrior (const parameters *p, int nP, likelihood *L) {
  double oldPrior=0, newPrior=0;
  int i;

  for(i=0; i<nP; i++) {
    double old=getParameter(p,i), new=getProposal(p,i);
    oldPrior-=1/expPar*old, newPrior-=1/expPar*new;
  }
  setPrior(L, oldPrior, newPrior);
/*   printPriorL(OUT, L); */

  return newPrior-oldPrior;
}


/*Reads a list of doubles*/
int readDoubles(FILE *fp, double *v, int buf) {
  int k=0;
  double d;
  char buff[100];
  while( (k<buf) && (fgets( buff, 100, fp ))) {
    /* If buff has more than just \n */
    if((strlen(buff)>1) && (sscanf(buff,"%lf", &d)))
      v[k++]=d;
  }
  return k;
}

/*Reads a dynamically allocated double matrix. Returns number of
  rows. */
int readMatrix(FILE *fp, double **A, int M, int N) {
  int k=0, l=0;
  double d;
  for(k=0; k<M; k++) {
    for(l=0; l<N; l++) {
      int err=fscanf(fp,"%lf", &d);
      if( err==EOF) {
	fprintf(ERR, "readMatrix(): File ended in line %i, before M=%i\n",
		k,  M);
	return k;
      }
      if(!err) {
	fprintf(ERR, "readMatrix(): Error reading line k=%i\n", k), exit(1);
      }
      A[k][l]=d;
    }
  }
  return k;
}

/* Returns true if data point is below threshold */
static int isOpen(double dat, double thresh) {
  return fabs(dat) > fabs(thresh);
}


/* Estimates the open sequences and the length of OCO sequences by
   thresholding and saves to NC and NO respectively - returns the
   length of the NC "histogram" */
int guessOClength(const double *data, int *NC, int *NO, int *n, double thresh, double *pC) {
  int i,j, maxClosed=0, open=0, closed=0;
  
  /*Get to first open state*/
  for(i=0; (i<(*n))  && !(isOpen(data[i], thresh)); i++);
  fprintf(OUT, "%i is first open state: %f\n",i, data[i]);

  /* Cut off if data does not end in open state */
  for(j=(*n)-1; j>=i && !(isOpen(data[j], thresh)); j--);
  (*n)=j;
  fprintf(OUT, "%i valid data points, data ends in: %f.\n",
	  (*n),data[(*n)-1]);

  while(i<(*n)) {
    /*count open states*/
    for(;(i<(*n)) && (isOpen(data[i],thresh)); i++,open++);
    /* fprintf(OUT, "%i: %i open states.\n", i, --open); */
    if( !(i<(*n) )) break;
    
    /*count closed states*/
    for(;(i<(*n)) && !(isOpen(data[i],thresh)); i++, closed++);
    if(closed>maxClosed) maxClosed=closed;
    /* fprintf(OUT, "%i closed events.\n", closed); */
    if(closed >0)
      NC[--closed]++, closed=0;
    else
      { fprintf(ERR, "guessOClength(): 0 closed?!?\n");break;}
  }
  /* fprintf(OUT, "guessOClength(): Finished\n"); */
  (*NO)=open;
  *pC=1-(double)open/(*n);
  fprintf(OUT, "guessOClength(): Open probability pO = %f\n",
	  1-(*pC));
  return maxClosed;
}

/* Converts A to NO and NC histogram: Po is in column 0 of A, length */
/*    of segment is in column 2. Returns maximum observed length of */
/*    consecutive closed events. Sets pC to estimated closed */
/*    probability. */
int modes2OClength(const double **A, int *NC, int *NO, int mA, double Pthresh, double *pC) {
  int nOpen=0;
  int nData=0;
  int i=0, maxClosed=0;
  int isO=A[i][0]>Pthresh;
  
  /* Check if first segment is 'open' */
  while(i<mA && !isO) i++, isO=A[i][0]>Pthresh;

  for(; i<mA; i++) {
    isO=A[i][0]>Pthresh;
    if(isO) {
      int open=(int)round(A[i][2]);
      fprintf(ERR, "modes2OClength(): %i: %i open events\n", i, open);
      nOpen+=open, nData+=open;
    }
    else {
      int closed=0;
      while( (i< mA) && (! isO) ) {
	closed+=(int)round(A[i][2]);
	fprintf(ERR, "modes2OClength(): %i: %i closed events\n", i, closed);
	i++;
	isO=A[i][0]>Pthresh;
      }
      /* Last segment must be thrown away because it does not end in
	 open event */
      if( i>=mA) return maxClosed;
      else { /* We know that isO is false, so the current segment
		contains open events */
	int open=(int)round(A[i][2]);
	fprintf(ERR, "modes2OClength(): %i: %i open events\n", i, open);
	nOpen+=open, nData+=open;
      }
      if(closed > maxClosed) maxClosed=closed;
      NC[--closed]++, nData+=closed;

    }
  }

  (*NO)=nOpen;
  *pC=1-(double)nOpen/nData;
  fprintf(OUT, "guessOClength(): nData=%i, pO = %f, pC=%f\n",
	  nData, 1-(*pC), *pC);
  return maxClosed;
}

/* Starting from column */
int* mp_matrix2OClength(FILE* datafp, int *nRows, int nCol,
			int *NO, int *maxClosed,
			double Pthresh, double *pC) {
  double **A=malloc((*nRows)*sizeof(*A));
  int i,j;
  int *NC;
  
  for(i=0; i<(*nRows); i++) A[i]=malloc(nCol*sizeof(double));

  fprintf(OUT, "mp_matrix2OClength(): Allocated memory for modesData\n");

  *nRows=readMatrix(datafp, A, *nRows, nCol);
  for(i=0; i<(*nRows); i++) {
    for(j=0; j<nCol; j++) {
      fprintf(OUT, "%g\t", A[i][j]);
    }
    fprintf(OUT,"\n");
  }
  fprintf(OUT, "mp_matrix2OClength(): Successfully read %i lines of modesData\n", (*nRows));
  /* Convert to nOpen/OCO*/
  /* Converts A to NO and NC histogram: Po is in column 0 of A, 
     length of segment is in column 2. */
  NC=calloc(NDATA, sizeof(*NC));
  (*maxClosed)=modes2OClength((const double **)A,NC, NO, (*nRows), Pthresh, pC);
  for(i=0; i<(*maxClosed); i++) {
    if(NC[i])	fprintf(OUT, "OC^nO: n=%i: %i times \n", i+1, NC[i]);
  }
  /*Free A*/
  for(i=0; i<(*nRows); i++) free(A[i]);
  free(A);

  return NC;
}


/*Generate O/C trace from thresholded open probabilities: Po is in column 0, length of segment is in column 2 */
void modes2Events(int *events, const double **A, double Pthreshold,
		  int nA) {

  int i, k=0, nOpen=0, nClosed=0, N;
  for(i=0; i<nA; i++) {
    /* type of event */
    int isO=(A[i][0]>Pthreshold);
    int L=(int)round(A[i][2]);
    int j;
    /* fprintf(OUT, "modes2Events(): %f > %f = %i\n", A[i][0], Pthreshold, isO); */
    for(j=0; j<L; j++, k++) events[k]=isO;
    if(isO) nOpen += L;
    else nClosed +=L;
  }
  N=nOpen+nClosed;
  fprintf(OUT, "modes2Events(): %i events, Po=%f, Pc=%f\n",
	  N, (double) nOpen/N, (double)nClosed/N);
}

/* Read data from file fp, threshold values for obtaining a sequence
   of open and closed events. Calculate sequences of closed events NC
   and calculate number of open states NO. Calculate closed
   probability pC. Return maximum length of closed sequence.*/
int processDataOneOpen(FILE *fp, int *NC, int *NO, int *nData, double thresh, double *pC) {
  double * data=malloc(NDATA*sizeof(double));
  int maxClosed=0;
  /* int i, j, open=0, closed=0; */
  fprintf(OUT, "Processing data...\n");
  readDoubles(fp, data,NDATA);
  /* NC=calloc(NDATA, sizeof(int)); */

  maxClosed=guessOClength(data, NC, NO, nData, thresh, pC);
  
  fprintf(OUT, "%i events read. Maximal closed sequence %i, pC=%.3f\n",
	  (*nData), maxClosed, *pC);
  free(data);

  return maxClosed;
}


/* Reads data from a file datafp with nRows rows. Threshold must be
   provided in thresh.  Returns the NC "histogram" and saves the
   number of open events to NO. Maximum number of closed events saved
   to maxClosed, ACTUAL number of rows read passed to nRows, total
   size of data set saved in nData, estimated closed probability saved
   in pC.*/
int* columnsData2NOC(const char *datafn, double thresh, int *NO, int *maxClosed,
		     int* nRows, int *nData, double *pC) {
  FILE *datafp=fopen(datafn, "r");
  nw_tokline **data=malloc((*nRows)*sizeof(*data));
  double *trace;
  int nCol;
  int i;
  int *NC;
  
  *nRows=nw_data_fpToTokline(datafp, data);
  fclose(datafp);
  
  fprintf(ERR, "columnsDataNOC(): Successfully read %i lines\n", *nRows);
  
  nCol=nw_data_nCol((const nw_tokline *)data[0]);
  trace=malloc(((*nRows)-1)*(nCol-1)*sizeof(*trace));
  fprintf(ERR, "%s", "columnsDataNOC(): Allocated trace\n");

  (*nData)=nw_data_convertColumnsToDouble((const nw_tokline **)data, 1, nCol-1,
				       1,  (*nRows)-1, trace);
  /* for(i=0; i<nData; i++) fprintf(OUT, "%lf\n", trace[i]); */
  fprintf(ERR, "%i doubles in trace\n", *nData);
  NC=calloc(*nData, sizeof(int));
  (*maxClosed)=guessOClength((const double *)trace, NC, NO, nData, thresh, pC);
  free(trace);
  for(i=0; i<(*nRows); i++) free(data[i]);
  free(data);

  return NC;
}

/* Thresholds nData points in data by threshold. Save bitstring to
   events. */
static void thresholdData(int *events, const double *data, double threshold,
		   int nData) {
  int i, nOpen=0;
  double Po, Pc;

  for(i=0; i<nData; i++) {
    int isO=isOpen(data[i], threshold);
    /* fprintf(ERR, "thresholdData(): data[%i]=%f (%i)\n", i, data[i], */
    /* 	    isOpen(data[i],threshold)); */
    /* fprintf(ERR, "%i %5f %i\n", i+1, data[i], isO); */
        /* fprintf(ERR, "%i\n",isOpen(data[i],threshold)); */
    if(isO) nOpen++;
    events[i]=isO;
  }
  
  fprintf(OUT, "Data processed and class sequence computed,");
  fprintf(OUT,"openEvents= %i, closedEvents=%i.\n", nOpen, nData-nOpen);
  Po=(double)nOpen/nData;
  Pc=1-Po;
  fprintf(OUT,"Estimated probabilities: Po = %f %%, Pc = %f %%\n",
	  Po*100,Pc*100);
}

/* Reads data from a file datafp with nRows rows. Threshold must be
   provided in thresh.  Returns events, an integer array of zeroes of
   ones. The ACTUAL number of rows read passed to nRows, total size of
   data set saved in nData.*/
int* mp_columnsData2Events(const char *datafn, double thresh,
			   int* nRows, int *nData) {
  FILE *datafp=fopen(datafn, "r");
  nw_tokline **data=malloc((*nRows)*sizeof(*data));
  double *trace;
  int nCol;
  int i;
  int *events;
  
  *nRows=nw_data_fpToTokline(datafp, data);
  fclose(datafp);
  
  fprintf(ERR, "mp_columnsData2Events(): Successfully read %i lines\n", *nRows);
  
  nCol=nw_data_nCol((const nw_tokline *)data[0]);
  trace=malloc(((*nRows)-1)*(nCol-1)*sizeof(*trace));
  fprintf(ERR, "%s", "mp_columnsData2Events(): Allocated trace\n");

  (*nData)=nw_data_convertColumnsToDouble((const nw_tokline **)data, 1, nCol-1,
					  1,  (*nRows)-1, trace);
  /* for(i=0; i<nData; i++) fprintf(OUT, "%lf\n", trace[i]); */
  fprintf(ERR, "%i doubles in trace\n", *nData);
  events=calloc((*nData), sizeof(int));
  thresholdData(events, trace, thresh, *nData);
  free(trace);
  for(i=0; i<(*nRows); i++) free(data[i]);
  free(data);

  return events;
}


/* Read model from file modelfile in double vector model */
double * readModel(const char *modelfile, int nStates) {
  FILE *fp;
  int i;
  double * model;
  
  fp = fopen(modelfile, "r");
  model = malloc(nStates*nStates*sizeof(double));
  doubleVec(fp, model);

  for(i=0; i<nStates*nStates; i++) {
    fprintf(OUT, "%f\t", model[i]);
    if( (i+1)%nStates == 0 ) fprintf(OUT, "\n");
  }
  fprintf(OUT, "\n");

  fclose(fp);

  return model;
}

/* Get Arguments */
void mp_getArgsOneOpen(char **argv, int argc,char **datafn,
		       char **modelfile,
		       int *nStates,
		       int *nIter,
		       double *delta,
		       int *seed,
		       double *thresh,
		       char **ratesFn,
		       char **statFn,
		       char **likeFn) {
  int nArgs=5;
  char *endptr;
  const char *usage="mp_OneOpen data model nStates nIterations [delta] [seed] [prefix]";
  const char *defaultRatesFn="rates.dat";
  const char *defaultStatFn="statDist.dat";
  const char *defaultLikelihoodName="likelihood.dat";
  char *prefix="piggy";
  
  if(argc < nArgs) {
    fprintf(ERR, "Usage: %s\n", usage), exit(1);
  }

  /* datafn=argv[1]; */
  /* printf("Data file:%s\n", datafn); */
  (*datafn)=argv[1];
  printf("Data file:%s\n", (*datafn));

  (*modelfile) = argv[2];
  printf("Model file:%s\n", (*modelfile));

  (*nStates) = strtol(argv[3], &endptr, 10);
  printf("Number of states: %i\n", (*nStates));

  (*nIter) = strtol(argv[4], &endptr, 10);
  printf("Number of iterations: %i\n", (*nIter));

  /* Optional arguments */
  if(argc>=nArgs+1) {
    (*delta) = strtod(argv[nArgs], &endptr);
    printf("delta: Step size of sampler: %f\n", (*delta));
  }

  if(argc>=nArgs+2) {
    (*seed) = strtol(argv[nArgs+1], &endptr, 10);
    printf("Seed for random number generator: %i\n", (*seed));
  }

  if(argc>=nArgs+3) {
    prefix = argv[nArgs+2];
    printf("Prefix for output files: %s\n", prefix);    
  }
  
  fprintf(OUT, "Threshold: %f\n", *thresh);

  /* Add prefix to filenames */
  (*ratesFn)=malloc(sizeof(char)*(strlen(*ratesFn)+strlen(prefix))+1);
  strcpy(*ratesFn, prefix);
  strcat(*ratesFn, defaultRatesFn);

  (*likeFn)=malloc(sizeof(char)*(strlen(prefix)+strlen(defaultLikelihoodName))+1);
  (*statFn)=malloc(sizeof(char)*(strlen(*statFn)+strlen(prefix)+1));

  strcpy(*statFn, prefix);
  strcat(*statFn, defaultStatFn);

  strcpy(*likeFn, prefix);
  strcat(*likeFn, defaultLikelihoodName);
  fprintf(OUT, "Output written to:\n");
  fprintf(OUT, "%s\n",*ratesFn);
  fprintf(OUT, "%s\n",*statFn);
  fprintf(OUT, "%s\n",*likeFn);
}

/* Get Arguments for N open states. */
void mp_getArgsNOpen(char **argv, int argc,char **datafn,
		     char **modelfile,
		     int *nStates,
		     int *nOpen,
		     int *nIter,
		     double *delta,
		     int *seed,
		     double *thresh,
		     char **ratesFn,
		     char **statFn,
		     char **likeFn) {
  int nArgs=6;
  char *endptr;
  const char *usage="mp_NOpen data model nStates nOpen nIterations [delta] [seed] [prefix]";
  const char *defaultRatesFn="rates.dat";
  const char *defaultStatFn="statDist.dat";
  const char *defaultLikelihoodName="likelihood.dat";
  char *prefix="piggyN";
  
  if(argc < nArgs) {
    fprintf(ERR, "Usage: %s\n", usage), exit(1);
  }

  /* datafn=argv[1]; */
  /* printf("Data file:%s\n", datafn); */
  (*datafn)=argv[1];
  printf("Data file:%s\n", (*datafn));

  (*modelfile) = argv[2];
  printf("Model file:%s\n", (*modelfile));

  (*nStates) = strtol(argv[3], &endptr, 10);
  printf("Number of states: %i\n", (*nStates));

  (*nOpen) = strtol(argv[4], &endptr, 10);
  printf("Number of open states: %i\n", (*nOpen));

  (*nIter) = strtol(argv[5], &endptr, 10);
  printf("Number of iterations: %i\n", (*nIter));

  /* Optional arguments */
  if(argc>=nArgs+1) {
    (*delta) = strtod(argv[nArgs], &endptr);
    printf("delta: Step size of sampler: %f\n", (*delta));
  }

  if(argc>=nArgs+2) {
    (*seed) = strtol(argv[nArgs+1], &endptr, 10);
    printf("Seed for random number generator: %i\n", (*seed));
  }

  if(argc>=nArgs+3) {
    prefix = argv[nArgs+2];
    printf("Prefix for output files: %s\n", prefix);    
  }
  
  fprintf(OUT, "Threshold: %f\n", *thresh);

  /* Add prefix to filenames */
  (*ratesFn)=malloc(sizeof(char)*(strlen(*ratesFn)+strlen(prefix))+1);
  strcpy(*ratesFn, prefix);
  strcat(*ratesFn, defaultRatesFn);

  (*likeFn)=malloc(sizeof(char)*(strlen(prefix)+strlen(defaultLikelihoodName))+1);
  (*statFn)=malloc(sizeof(char)*(strlen(*statFn)+strlen(prefix)+1));

  strcpy(*statFn, prefix);
  strcat(*statFn, defaultStatFn);

  strcpy(*likeFn, prefix);
  strcat(*likeFn, defaultLikelihoodName);
  fprintf(OUT, "Output written to:\n");
  fprintf(OUT, "%s\n",*ratesFn);
  fprintf(OUT, "%s\n",*statFn);
  fprintf(OUT, "%s\n",*likeFn);
}


/* count positive entries in Q */
int countPositive(const gsl_matrix *Q) {
  int i, j, nPos=0;

  for(i=0; i<Q->size1; i++) {
    for(j=0; j<Q->size2; j++) {
      if(gsl_matrix_get(Q, i, j) > 0) nPos++;
    }
  }

  return nPos;
}


/* converts Markov model to parameters: qIJ in first half, qJI in
   second half, starting from index i0 */
/* Q cannot be const because org[k] pointers are not const */
void matrix2ParameterSubdiag(gsl_matrix *Q, gsl_matrix *Qp, 
			     parameters *p, int nPar,int *qInd) {
  int i, j, k;
  int centre =nPar/2;

  /* for(i=0, k=0; i<Q->size1; i++) { */
  /*   /\*subdiagonal!*\/ */
  /*   for(j=0; j<i; j++) { */
  for(j=0, k=0; j<Q->size2; j++) {
    /*subdiagonal!*/
    for(i=0; i<j; i++) {
    /* for(j=i; j<Q->size2 ; j++) { */
      if(gsl_matrix_get(Q, i, j) > 0) {
	char * name1 = malloc( (strlen("q_{ij}")+1)*sizeof(char));
	char * name2 = malloc( (strlen("q_{ij}")+1)*sizeof(char));

	sprintf(name1, "q_{%i%i}", i+1, j+1);
	setParameterName(p, k, name1);
	setParameterPtr(p, k, gsl_matrix_ptr(Q, i, j));
	setProposalPtr(p, k, gsl_matrix_ptr(Qp, i, j));

	sprintf(name2, "q_{%i%i}", j+1, i+1);
	setParameterName(p, centre + k, name2);
	setParameterPtr(p, centre+k, gsl_matrix_ptr(Q, j, i));
	setProposalPtr(p, centre + k, gsl_matrix_ptr(Qp, j, i));

	qInd[k]=i, qInd[k+centre]=j;
	k++;

	/* free strings because setParameterName copies anyway*/
	free(name1), free(name2);
      }
    }
  }	
}

/* converts Markov model to parameters: qIJ is saved at even, qJI at odd
   indices starting from k0 */
void matrix2ParameterPairsK0(gsl_matrix *Q, gsl_matrix *Qp, 
			     parameters *p, 
			     int k0, int nPar,int *qInd) {
  int i, j, k;

  for(i=0, k=k0; i<Q->size1; i++) {
    /*subdiagonal!*/
    /* for(j=0; j<i; j++) { */
    for(j=i; j<Q->size2; j++) {
      if(gsl_matrix_get(Q, i, j) > 0) {
	char * name1 = malloc( (strlen("q_{ij}")+1)*sizeof(char));
	char * name2 = malloc( (strlen("q_{ij}")+1)*sizeof(char));

	sprintf(name1, "q_{%i%i}", i+1, j+1);
	setParameterName(p, k, name1);
	setParameterPtr(p, k, gsl_matrix_ptr(Q, i, j));
	setProposalPtr(p, k, gsl_matrix_ptr(Qp, i, j));
	qInd[k]=i;
	k++;

	sprintf(name2, "q_{%i%i}", j+1, i+1);
	setParameterName(p, k, name2);
	setParameterPtr(p, k, gsl_matrix_ptr(Q, j, i));
	setProposalPtr(p, k, gsl_matrix_ptr(Qp, j, i));

	qInd[k]=j;
	k++;

	/* free strings because setParameterName copies anyway*/
	free(name1), free(name2);
      }
    }
  }	
}

/* converts Markov model to parameters: qIJ at even, qJI in odd
   indices. Q cannot be const because org[k] pointers are not const */ 
void matrix2ParameterPairs(gsl_matrix *Q, gsl_matrix *Qp, 
			   parameters *p, 
			   int nPar,int *qInd) {
  matrix2ParameterPairsK0(Q, Qp, p, 0, nPar, qInd);
}

/*output nRates rate constants starting from index k0 and additionally
  calculate stationary distribution*/
void QandPhioutK0(FILE *ratesfp, int iteration, 
		  const parameters *rates, int k0, int nP, 
		  FILE *statfp, gsl_matrix *Q) {
  gsl_vector *phi=gsl_vector_alloc(Q->size1);

  makeConservative(Q);
  statMarkovVector(Q, phi);

  fprintf(statfp, "%i\t", iteration);
  vprintfmt(statfp, phi,"%g");
  fflush(statfp);

  generalOutputK0(ratesfp, iteration, rates, k0, nP);  

  gsl_vector_free(phi);
}

/*output rates and additionally calculate stationary distribution*/
void QandPhiout(FILE *ratesfp, int iteration, 
		const parameters *rates, int nP, 
		FILE *statfp, gsl_matrix *Q) {
  QandPhioutK0(ratesfp,  iteration, rates,  0,  nP, statfp, Q);
}


/* scale second half of parameters using the stationary
   distribution. */
static void scaleRates(parameters *rates, int nRates, 
		       const double *phi, const int *qInd) {
  int i, centre=nRates/2;

  for(i=0; i<centre; i++) {
    double qIJ = getProposal(rates, i), qJI;
    double phiI=phi[qInd[i]], phiJ=phi[qInd[i+centre]];

/*     fprintf(OUT, "phi[%i] = %f, phi[%i] = %f, qIJ = %f, qJI = %f\n", */
/*     	    qInd[i], phiI, qInd[i+centre], phiJ, qIJ,  getProposal(rates, i+centre)); */

    if(phiJ < EPS) {
      fprintf(ERR, "scaleRates(): Division by small phi(%i) = %g. Exiting...\n", 	      qInd[i+centre],phiJ);
      exit(1);
    }
    qJI=phiI/phiJ*qIJ;
    /* fprintf(OUT, "new qJI = %f\n", qJI); */
    setProposal(rates, i+centre, qJI);
  }
}

/* Checks if v contains a small entry */
static int smallEntry(const double *v, int n){
  int i;
  for (i=0; i<n; i++) {
    if(v[i] < EPS) return 1;
  }
  return 0;
}

/* Gibbs sampling strategy for sampling rates and rate constants */
/* alpha is the set of hyperparameters for the rate constants */
/* ONE COULD GIVE a vector which is a proposal for the stationary
   distribution as an argument */
void samplePhiandQ(const gsl_rng *r, 
		   parameters *p, int nPar,
		   const double *alpha,
		   gsl_matrix *Q,
		   const int *qIndices) {

  int nStates=Q->size1;
  int nRates=nPar/* -nStates */;
  /* gsl_vector *statDist=gsl_vector_alloc(nStates); */
  double *statDist=malloc(nStates*sizeof(double));

  int count, i;
  const int maxCount=10000;

  /* fprintf(OUT, "proposePhi(): Change phi...\n"); */
  count =1;
  do {
    gsl_ran_dirichlet (r, nStates, alpha, statDist);
  }  /* while( (fabs(gsl_vector_min(phiV)) < EPS) && (++count <= maxCount)); */
  while( (smallEntry(statDist, nStates)) && (++count <= maxCount));
  if (count >= maxCount) {
    fprintf(ERR, "proposePhi(): no positive sample of stationary distribution\nafter %i iterations. Exiting...\n", maxCount);
    exit(1);
  }

  /*Set parameter to stattionary distribution*/
  for(i=0; i<nStates; i++) {
    setProposal(p, nRates+i, statDist[i]);
  }

/*   if(! gsl_vector_ispos(&phiPView.vector)) { */
/*     fprintf(ERR, "Sampled non-positive vector. Exiting...\n"); */
/*     vprint(ERR, &phiPView.vector); */
/*     exit(1); */
/*   } */

/*   fprintf(OUT, "proposePhi(): phi = "); */
/*   for(i=0; i<nStates; i++) { */
/*       fprintf(OUT, "%f\t", phi[i]); */
/*     } */
/*   fprintf(OUT, "\n"); */

  /* Sample first half of rate constants --> parP*/
/*   sampleConstrained(r, p, nRates/2,NULL, 0); */
  proposeConstrainedAB(r, p, 0, nRates);

/*   fprintf(OUT, "samplePhiandQ(): Parameter:\n"); */
/*   printParameters(OUT, p,nRates+nStates); */
/*   fprintf(OUT, "samplePhiandQ(): Proposal:\n"); */
/*   printProposal(OUT, p,nRates+nStates); */

  /* scale remaining half using the stationary distribution */
  scaleRates(p, nRates, statDist, qIndices);

  free(statDist);

/*   fprintf(OUT, "samplePhiandQ(): Parameter:\n"); */
/*   printParameters(OUT, p,nRates+nStates); */
/*   fprintf(OUT, "samplePhiandQ(): Proposal:\n"); */
/*   printProposal(OUT, p,nRates+nStates); */

/*   exit(1); */
}

/* Sample steps parameter according to constraints starting from n0 */
void sampleN(const gsl_rng *r, 
	     parameters *p, int n0, int steps) {
  proposeConstrainedAB(r, p, n0, steps);
}

/* Checks if matrix is symmetrix */
int symmetricQ(gsl_matrix* A) {
  int i, j;
  for(i=0; i<A->size1; i++) {
    for(j=i+1;j<A->size2; j++) {
      if(fabs(gsl_matrix_get(A,i,j)-gsl_matrix_get(A,j,i)) > EPS) return 0;
    }
  }
  return 1;
}


/* Checks if a model Q satisfies the detailed balance condition */
int checkDetailedBalance(const gsl_matrix *Q) {
  gsl_matrix *pi=gsl_matrix_alloc(Q->size1,Q->size2);
  gsl_matrix *D=gsl_matrix_calloc(Q->size1,Q->size2);
  gsl_matrix *S=gsl_matrix_alloc(Q->size1,Q->size2);
  int i;
  int db;
  const int showSym=1;

  printMat(OUT, Q);
  pi=statMarkov(Q);

  for(i=0; i<pi->size2; i++)
    gsl_matrix_set(D, i, i, gsl_matrix_get(pi, 0,i));
    
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1, D, Q, 0, S);

  if(showSym) {
    fprintf(OUT, "Reparametrisation of Q with stationary distribution:\n");
    printMat(OUT, S);
  }
  fprintf(OUT, "Detailed balance? ");
  db=symmetricQ(S);
  if(db) fprintf(OUT,"TRUE\n");
  else fprintf(OUT,"FALSE\n");
  gsl_matrix_free(S), gsl_matrix_free(D), gsl_matrix_free(pi);

  return db;
}


/* Enforces detailed balance for a given matrix with cycle */
void forceDetailedBalance(gsl_matrix *Qp,
			  const int *cycle, int nCycle) {
  int i,j;
  double numerator=1/* gsl_matrix_get(Qp,cycle[nCycle-1], cycle[0]); */;
  double denominator=1;
  double qIJ,qJI;

  for(i=0, j=1; j< nCycle; i++, j++) {
    qIJ=gsl_matrix_get(Qp, cycle[i], cycle[j]);
    qJI=gsl_matrix_get(Qp, cycle[j], cycle[i]);
    numerator *= qIJ, denominator *= qJI;    
  }
  if(denominator < EPS) {
    fprintf(ERR,"sampleQdetailedBalance: Division by %g - too small.\n", denominator);
    exit(1);
  }
  qIJ=gsl_matrix_get(Qp,cycle[nCycle-1], cycle[0]);
  qJI=gsl_matrix_get(Qp,cycle[0], cycle[nCycle-1]);
  if(numerator*qIJ<denominator) {
    gsl_matrix_set(Qp, cycle[0], cycle[nCycle-1], qIJ*numerator/denominator);
  }
  else {
    gsl_matrix_set(Qp, cycle[nCycle-1], cycle[0], qJI*denominator/numerator);
  }
}

/* Array containing cycle indices in the right sequence*/
/* int *cycle, nCycle; */

/* Samples rate constants and enforces detailed balance for a cycle,
   given by an array of indices*/
/* void sampleQdetailedBalance(const gsl_rng *r,  */
/* 			    parameters* p, int nPar, */
/* 			    intparameters *ip, int nIPar) { */
/*   sampleConstrained(r, p, nPar, ip, nIPar); */
/*   forceDetailedBalance(Qp, cycle, nCycle); */
/* } */
