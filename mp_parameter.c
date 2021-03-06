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

#include <stdio.h>
#include <string.h>
#include <math.h>

#ifndef HPC
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#else
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#endif

#include "mp_parameter.h"

#define OUT stdout
#define ERR stderr

#define EPS 1e-12

#define FLIP(x) ( (x) != 0) ? 0 : 1




/* Dynamically allocates room for s plus additional n characters */
static char *strnsave(const char * s, size_t n) {
  char *p = malloc(strlen(s) + n + 1);
  if(!p) fputs("strsave: out of memory\n", ERR), exit(1);

  return strcpy(p, s);
}

/* Dynamically allocates room for s plus additional n characters */
static char *strsave(const char * s) {
  return strnsave(s,0);
}

/* gives i-th of the proposals */
double getParameter(const parameters *p, int i) {
  return *(p[i].org);
}

/* gives i-th of the proposal */
double getProposal(const parameters *p, int i) {
  return *(p[i].prop);
}

/* gives i-th sampling width of the parameters */
double getWidth(const parameters *p, int i) {
  /* parameters par=p[i]; */
  /* printFullParameters(OUT,&par,1); */
  return p[i].width;
}

/* set i-th sampling width of the parameters */
void setWidth(parameters *p, int i, double d) {
  p[i].width=d;
}

/* set sampling widths of the parameters */
void setWidths(parameters *p, const double *d, int nPar) {
  int i;
  for(i=0; i<nPar; i++) setWidth(p, i, d[i]);
}

/* gives i-th sampling width of the parameters */
double getMin(const parameters *p, int i) {
  return p[i].min;
}

/* set i-th sampling width of the parameters */
void setMin(parameters *p, int i, double d) {
  p[i].min=d;
}

/* gives i-th sampling width of the parameters */
double getMax(const parameters *p, int i) {
  return p[i].max;
}

/* set i-th sampling width of the parameters */
void setMax(parameters *p, int i, double d) {
  p[i].max=d;
}

/* set i-th parameter to d */
void setParameter(parameters *p, int i, double d) {
  *(p[i].org) = d;
}

/* set parameters to components of d */
void setParameters(parameters *p, int pn0,
		   double* const d, int dn0, int nD) {
  int i;
  for(i=0; i<nD; i++)
    setParameter(p, pn0+i, d[dn0+i]);
}

/* set i-th parameter to d */
void setParameterPtr(parameters *p, int i, double *d) {
  p[i].org = d;
}

/* set i-th proposal to d */
void setProposal(parameters *p, int i, double d) {
  *(p[i].prop)=d;
}

/* set proposals to components of d */
void setProposals(parameters *p, int pn0,
		  double* const d, int dn0, int nD) {
  int i;
  for(i=0; i<nD; i++)
    setProposal(p, pn0+i, d[dn0+i]);
}

/* set i-th proposal to d */
void setProposalPtr(parameters *p, int i, double *d) {
  p[i].prop=d;
}

/* copies values from original parameters to proposal */
void copyOrg2Prop(parameters *p, int nPar) {
  int i;
  for (i=0; i<nPar; i++) setProposal(p, i, getParameter(p,i));
}

/* copies values from proposal to original parameters */
void copyProp2OrgK0(parameters *p, int k0, int nPar) {
  int i;
  for (i=k0; i<k0+nPar; i++) {
    /* fprintf(OUT, "copyProp2Org(): p[%i]=%f, proposal=%f\n", */
    /* 	    i,getParameter(p,i), getProposal(p,i)); */
    setParameter(p, i, getProposal(p,i));
  }
}

/* copies values from proposal to original parameters */
void copyProp2Org(parameters *p, int nPar) {
  copyProp2OrgK0(p, 0,  nPar);
}

/* Generic "swapper" for pointers */
static void swap(void **x, void **y) {
	void *t = *x;
	*x = *y;
	*y = t;
}

static void swapPointers(void * pOne, void * pTwo)
{
  void * pTemp;
  pTemp = pOne;
  pOne = pTwo;
  pTwo = pTemp;
}

/* Swap parameters i and j */
void swapParameters(parameters *p, int i, int j) {
  /*swap name */
  swap((void**)&p[i].name,(void**)&p[j].name);
  /*swap org */
  swap((void**)&p[i].org,(void**)&p[j].org);
  /*swap prop */
  swap((void**)&p[i].prop,(void**)&p[j].prop);

  /* swap width */
  swapPointers(&p[i].width,&p[j].width);
  /* swap width */
  swapPointers(&p[i].min,&p[j].min);
  /* swap width */
  swapPointers(&p[i].max,&p[j].max);
}

/* initialise parameters with random numbers, uniformly distributed
   over [0,scale] */
void initUniformScale(const gsl_rng *r, parameters *par, 
		      double scale, int nPar) {
  int i;
  for(i=0; i<nPar; i++) setParameter(par, i, scale*gsl_rng_uniform(r));
  copyOrg2Prop(par, nPar);
}


/* initialise parameters with random numbers, uniformly distributed
   over [0,1] */
void initUniform(const gsl_rng *r, parameters *par, int nPar) {

  initUniformScale(r, par, 1, nPar);
}

static int cmpDouble(const void *d1, const void *d2) {
  return ((*(double*)d1)-(*(double*)d2))>0;
}


/* initialise parameters with random numbers, uniformly distributed
   over {i0,..., n-1}, sort in ascending order */
void initUniformAscending(const gsl_rng *r, parameters *par, int nPar) {
  int i;
  double *sortedPars=malloc(nPar*sizeof(double));

  for(i=0; i<nPar; i++) {
    sortedPars[i]=gsl_rng_uniform(r);
  }
  qsort(sortedPars, nPar, sizeof(double), cmpDouble);
  for(i=0; i<nPar; i++) {
    setParameter(par, i, sortedPars[i]);
  }

  copyOrg2Prop(par, nPar);
}



/*initialises parameter with default values for width, min and max */
void initialiseParameters(parameters *p, 
			  double width, double min, double max, int nPar) {
  int i;

  for(i=0; i<nPar; i++) {
    char *name = malloc(strlen("p_{00}")+1);
    sprintf(name,"p_{%i}", i);
    setParameterName(p, i, name);
    /* p[i].org=malloc(sizeof(double)); */
    /* p[i].prop=malloc(sizeof(double)); */
    setWidth(p, i, width);
    setMin(p, i, min);
    setMax(p, i, max);
    /* free string because setParameterName copies anyway*/
    free(name);
  }
}

/* Copies parameter new to dest, reserves memory for org/prop when
   necessary */
void parameterMemcpy(parameters *dest, const parameters *new, int nPar){
  int i;
  for(i=0; i<nPar; i++) {
    dest[i].name=new[i].name;
    /* if(!dest[i].org) dest[i].org=malloc(sizeof(double)); */
    dest[i].org=new[i].org;
    /* setParameter(dest,i, getParameter(new,i)); */
    /* if(!dest[i].prop) dest[i].prop=malloc(sizeof(double)); */
    dest[i].prop=new[i].prop;
    /* setProposal(dest,i, getProposal(new,i)); */
    setWidth(dest,i,getWidth(new,i));
    setMin(dest,i,getMin(new,i));
    setMax(dest,i,getMax(new,i));
  }
}

/* returns name of i-th parameter */
char * getParameterName(const parameters *p, int i) {
  return p[i].name;
}

/* sets name of i-th parameter */
void setParameterName(parameters *p, int i, char * const name) {
  p[i].name=strsave(name);
}

/* sets names of parameters, starting index at k0 */
void setParameterNamesK0(parameters *p, char * const prefix, int k0, int n) {
  const int digits=3;
  int i;
  for(i=k0; i<k0+n; i++) {
    p[i].name=malloc((strlen(prefix)+digits+1)*sizeof(char));
    sprintf(p[i].name, "%s%03i", prefix, i);
  }
}

/* sets names of parameters */
void setParameterNames(parameters *p, char * const prefix, int n) {
  setParameterNamesK0(p,prefix, 0, n);
}

/*prints parameter names to file fp */
void printParameterNames(FILE *fp, const parameters *p, int nPar) {
  int i;

  for(i=0; i<nPar-1; i++) {
    fprintf(fp, "%s\t", getParameterName(p, i));
  }
  fprintf(fp, "%s\n", getParameterName(p, nPar-1));
}

/* saves parameters in file fp, starting from index k0 */
void printParametersK0(FILE *fp, const parameters *p, int k0, int nPar) {
  int i;
  for(i=0; i<nPar-1; i++) {
    fprintf(fp, "%f\t", getParameter(p, k0+i));
  }
  fprintf(fp, "%f\n", getParameter(p,k0+nPar-1));
}

/* saves parameters in file fp */
void printParameters(FILE *fp, const parameters *p, int nPar) {
  printParametersK0(fp, p, 0, nPar);
}

/* saves parameters in file fp with all additional information */
void printFullParameters(FILE *fp, const parameters *p, int nPar) {
  int i;

  for(i=0; i<nPar; i++) {
    fprintf(fp, "%s=%f (%f), w=%f, min=%f, max=%f\n", 
	    getParameterName(p, i), getParameter(p, i), getProposal(p, i),
	    p[i].width, p[i].min, p[i].max);
  }
}

/* saves proposal in file fp */
void printProposal(FILE *fp, const parameters *p, int nPar) {
  int i;

  for(i=0; i<nPar-1; i++) {
    fprintf(fp, "%f\t", getProposal(p, i));
  }
  fprintf(fp, "%f\n", getProposal(p,nPar-1));    
}

/* Initialises Parameters with pointers to elements in double array */
void initialisePorg(parameters *p, double *d, int n) {
  int i;
  for(i=0; i<n; i++) {
    p[i].org=d+i;
  }
}

/* Initialises parameter proposals with pointers to elements in double array */
void initialisePprop(parameters *p, double *d, int n) {
  int i;
  for(i=0; i<n; i++) {
    p[i].prop=d+i;
  }
}

#ifndef PROP
/*Returns next value for uniform random walk*/
static double walk(const gsl_rng *r, double width) {
  return (gsl_rng_uniform(r)-0.5)*width;
}

/* Walk until positive */
static double walkPositive(const gsl_rng *r, 
			   const double par, const double delta) {
  double w;
  do {
    w = walk(r, delta);
  } while(par + w <= 0);
  return w;
}

static void changePositive(const gsl_rng *r, parameters *p, 
			   int i, double width) {
  double par=getParameter(p, i);
  double w=walkPositive(r, par, width);
  setProposal(p, i, par+w);
}

/* Propose parameters */
void proposePositiveAB(const gsl_rng *r, parameters *p, int i0, int steps) {
  int i, n=i0+steps;

  for(i=i0; i<n; i++) {
    changePositive(r, p, i,  getWidth(p,i));
  }
}

/* check if proposal is in the right range */
int isProposalMinMax(const parameters *p, int nPar) {
  int i;

  for(i=0; i<nPar; i++) {
    double prop = getProposal(p, i);
    if( (prop < getMin(p,i)) || (prop > getMax(p,i)) ) {
/*       fprintf(ERR, "negative proposal found:\n"); */
      return 0;
    }
  }
  return 1;
}


static void changeConstrained(const gsl_rng *r, parameters *p, 
			      int i) {
  double par=getParameter(p, i);
  double w;
  double min = getMin(p,i), max = getMax(p,i);
  do {
    double width = getWidth(p,i);
    w = walk(r, width);
    /* fprintf(OUT, "changeConstrained(): i = %i, width = %f, par = %f, w = %f, min=%f, max=%f\n", */
    /* 	    i, width, par, w, min, max); */
  } while( (par + w < min) || (par + w > max) );
  setProposal(p, i, par+w);
  /* fprintf(OUT, "Proposal=%f\n",getProposal(p,i)); */
}

/* Propose parameters */
void proposeConstrainedAB(const gsl_rng *r, parameters *p, int i0, int steps) {
  int i;
  for(i=i0; i<i0+steps; i++) {
    changeConstrained(r, p, i);
  }
  /* fprintf(OUT, "proposeConstrainedAB():\n"); */
  /* printParameters(OUT, p, steps); */
  /* fprintf(OUT, "proposeConstrainedAB():\n"); */
  /* printProposal(OUT, p, steps); */
}

/* Propose parameters: Ratio of new and old is log(uNew/uOld)=u, where
   u is uniformly distributed over (-0.5, 0.5). */
void proposeExpAB(const gsl_rng *r, parameters *p, int i0, int steps) {
  int i;
  for(i=i0; i<i0+steps; i++) {
    double u=gsl_rng_uniform(r)-0.5;
    double oldP=getParameter(p,i);
    /* exp of uniform random number in interval ]-0.5,0.5 [*/
    /* fprintf(OUT, "proposeExpAB(): oldP = %f, newP=%f\n", oldP, oldP*exp(u)); */
    setProposal(p, i, oldP*exp(u));
    /* changeConstrained(r, p, i); */
  }
}



/* Does proposal match constraints ? */
int matchConstraintsK0(const parameters *p, int k0, int nPar) {
  int i;
  for(i=k0; i<k0+nPar; i++) {
    double prop = getProposal(p, i);
    double min = getMin(p, i);
    double max = getMax(p, i);
    if( (prop < min) || (prop > max)) {
      return 0;
    }
  }
  return 1;
}

/* Does proposal match constraints ? */
int matchConstraints(const parameters *p, int nPar) {
  return matchConstraintsK0(p, 0, nPar);
}

/* check if proposal is positive */
int isProposalPositive(const parameters *p, int nPar) {
  int i;

  for(i=0; i<nPar; i++) {
    double prop = getProposal(p, i);
    if( prop < EPS) {
/*       fprintf(ERR, "negative proposal found:\n"); */
      return 0;
    }
  }
  return 1;
}
#endif
/******************* End: double parameters ********************/

/****** int parameters *************/
/* gives i-th of the proposals */
int getIntParameter(const intparameters *p, int i) {
  return *(p[i].org);
}

/* gives i-th of the proposal */
int getIntProposal(const intparameters *p, int i) {
  return *(p[i].prop);
}

/* gives i-th of the parameters */
int getIntWidth(const intparameters *p, int i) {
  return p[i].width;
}

/* set i-th sampling width of the parameters */
void setIntWidth(intparameters *p, int i, int d) {
  p[i].width=d;
}

/* set sampling widths of the parameters */
void setIntWidths(intparameters *p, const int *d, int nPar) {
  int i;
  for(i=0; i<nPar; i++) setIntWidth(p, i, d[i]);
}

/* gives i-th sampling width of the parameters */
int getIntMin(const intparameters *p, int i) {
  return p[i].min;
}

/* set i-th sampling width of the parameters */
void setIntMin(intparameters *p, int i, int d) {
  p[i].min=d;
}

/* gives i-th sampling width of the parameters */
int getIntMax(const intparameters *p, int i) {
  return p[i].max;
}

/* set i-th sampling width of the parameters */
void setIntMax(intparameters *p, int i, int d) {
  p[i].max=d;
}

/*initialises parameter with default values for width, min and max */
void initialiseIntParameters(intparameters *p, int width, int min, int max, int nPar) {
  int i;

  for(i=0; i<nPar; i++) {
    /* p[i].org=malloc(sizeof(int)); */
    /* p[i].prop=malloc(sizeof(int)); */
    setIntWidth(p, i, width);
    setIntMin(p, i, min);
    setIntMax(p, i, max);
  }
}

/* set i-th parameter to d */
void setIntParameter(intparameters *p, int i, int d) {
  *(p[i].org)=d;
}

/* set i-th proposal to d */
void setIntParameterPtr(intparameters *p, int i, int *d) {
  p[i].org=d;
}

/* set parameters to components of d */
void setIntParameters(intparameters *p, const int* d, int nD) {
  int i;
  for(i=0; i<nD; i++)
    setIntParameter(p, i, d[i]);
}

/* set i-th proposal to d */
void setIntProposal(intparameters *p, int i, int d) {
  *(p[i].prop)= d;
}

/* set i-th proposal to d */
void setIntProposalPtr(intparameters *p, int i, int *d) {
  p[i].prop=d;
}

/* copies values from original parameters to proposal */
void copyIntOrg2Prop(intparameters *p, int nPar) {
  int i;
  for (i=0; i<nPar; i++) setIntProposal(p, i, getIntParameter(p,i));
}

/* resets all proposals to original parameter values */
void resetProposals(parameters *p, int nPar,
		    intparameters *ip, int nIPar) {
  copyOrg2Prop(p, nPar);
  copyIntOrg2Prop(ip, nIPar);
}

/* copies values from proposal to original parameters */
void copyIntProp2Org(intparameters *p, int nPar) {
  int i;
  for (i=0; i<nPar; i++) setIntParameter(p, i, getIntProposal(p,i));
}

/* Returns uniformly distributed integers from {i0, ..., i0+n-1} */
/* static */ int intUniform(const gsl_rng *r, int i0, int n) {
  int w=n;
  int res;

  if(w <=0) {
    fprintf(ERR, "intUniform(): Not possible to sample between %i and %i\n",
	    i0, n);
    exit(1);
  }
  res=i0+gsl_rng_uniform_int (r, w);

  return res;
}

/* initialise parameters with random numbers, uniformly distributed
   over {i0,..., n-1} */
void initIntUniform(const gsl_rng *r, intparameters *par, int i0, int n, int nPar) {
  int i;
  for(i=0; i<nPar; i++) {
    setIntParameter(par, i, intUniform(r, i0, n));
  }

  copyIntOrg2Prop(par, nPar);
}

static int cmpInt(const void *i1, const void *i2) {
  return (*(int*)i1)-(*(int*)i2);
}


/* proposal, uniformly distributed over {i0,..., n-1}, sort in
   ascending order */
void proposeIntUniformAscending(const gsl_rng *r, intparameters *par, int i0, int n, int nPar) {
  int i, *sortedPars=malloc(nPar*sizeof(int));

  /* fprintf(OUT,"proposeIntUniformAscending():\n"); */
  for(i=0; i<nPar; i++) {
    sortedPars[i]=intUniform(r, i0, n);
  }
  qsort(sortedPars, nPar, sizeof(int), cmpInt);
  /* fprintf(OUT,"Sorting completed...\n"); */
  for(i=0; i<nPar; i++) {
    setIntProposal(par, i, sortedPars[i]);
  }
  /* fprintf(OUT,"Proposals set...\n"); */
  free(sortedPars);
}

/* initialise parameters with random numbers, uniformly distributed
   over {i0,..., n-1}, sort in ascending order */
void initIntUniformAscending(const gsl_rng *r, intparameters *par, int i0, int n, int nPar) {
  int i, *sortedPars=malloc(nPar*sizeof(int));

  for(i=0; i<nPar; i++) {
    sortedPars[i]=intUniform(r, i0, n);
  }
  qsort(sortedPars, nPar, sizeof(int), cmpInt);
  for(i=0; i<nPar; i++) {
    setIntParameter(par, i, sortedPars[i]);
  }

  copyIntOrg2Prop(par, nPar);
}

/* Copies paramerter new to dest */
void intParameterMemcpy(intparameters *dest, const intparameters *new, int nPar){
  int i;

  for(i=0; i<nPar; i++) {
    dest[i].name=new[i].name;

    /* if(!dest[i].org) dest[i].org=malloc(sizeof(int)); */
    dest[i].org=new[i].org;

    /* if(!dest[i].prop) dest[i].prop=malloc(sizeof(int)); */
    dest[i].prop=new[i].prop;
    setIntWidth(dest,i,getIntWidth(new,i));
    setIntMin(dest,i,getIntMin(new,i));
    setIntMax(dest,i,getIntMax(new,i));
  }
}


/* returns name of i-th parameter */
char * getIntParameterName(const intparameters *p, int i) {
  return p[i].name;
}

/* returns name of i-th parameter */
void setIntParameterName(intparameters *p, int i, char *name) {
  p[i].name= strsave(name);
}

/* sets names of parameters */
void setIntParameterNames(intparameters *ip, char * const prefix, int n) {
  int digits=2;
  int i;
  for(i=0; i<n; i++) {
    ip[i].name=malloc((strlen(prefix)+digits+1)*sizeof(char));
    sprintf(ip[i].name, "%s%03i", prefix, i);
  }
}


/*prints parameter names to file fp */
void printIntParameterNames(FILE *fp, const intparameters *p, int nPar) {
  int i;

  for(i=0; i<nPar-1; i++) {
    fprintf(fp, "%s\t", getIntParameterName(p, i));
  }
  fprintf(fp, "%s\n", getIntParameterName(p, nPar-1));
}

/* saves parameters in file fp */
void printIntParameters(FILE *fp, const intparameters *p, int nPar) {
  int i;

  for(i=0; i<nPar-1; i++) {
    fprintf(fp, "%i\t", getIntParameter(p, i));
  }
  fprintf(fp, "%i\n", getIntParameter(p,nPar-1));
    
}

/* Prints constraints for each int parameter: width of random walk,
   min, max */
void printIntConstraints(FILE *fp, const intparameters *p, int nPar) {
  int i;
  for(i=0; i<nPar; i++) {
    fprintf(OUT, "%i: width=%i, min=%i, max=%i\n", i,
  	    getIntWidth(p, i), getIntMin(p,i), getIntMax(p, i));
  }
}

/* saves proposal in file fp */
void printIntProposal(FILE *fp, const intparameters *p, int nPar) {
  int i;

  for(i=0; i<nPar-1; i++) {
    fprintf(fp, "%i\t", getIntProposal(p, i));
  }
  fprintf(fp, "%i\n",getIntProposal(p, nPar-1));
    
}

/* saves parameters in file fp with all additional information */
void printIntFullParameters(FILE *fp, const intparameters *ip, int nIPar) {
  int i;

  for(i=0; i<nIPar; i++) {
    fprintf(fp, "%s=%i (%i), w=%i, min=%f, max=%f\n", 
 	    getIntParameterName(ip, i), getIntParameter(ip, i), getIntProposal(ip, i),
	    ip[i].width, ip[i].min, ip[i].max);
  }
}

/* Initialises Parameters with pointers to elements in double array */
void initialiseIntPorg(intparameters *ip, int *k, int n) {
  int i;
  for(i=0; i<n; i++) {
    ip[i].org=k+i;
  }
}

/* Initialises parameter proposals with pointers to elements in double array */
void initialiseIntPprop(intparameters *ip, int *k, int n) {
  int i;
  for(i=0; i<n; i++) {
    ip[i].prop=k+i;
  }
}


/* /\* Generates either sample of -width or +width *\/ */
/* static intUpDown(const gsl_rng *r, int width) { */
/*   int w = intUniform(r, 0,width); */

/*   if(gsl_rng_uniform(r) < 0.5) return w; */

/*   return -w; */
/* } */

/*Returns next value for uniform random walk*/
static double intSymmWalk(const gsl_rng *r, int width) {
  return intUniform(r, -width, 2*width+1);
}

/* Walk until positive */
static double intSymmWalkPositive(const gsl_rng *r, 
				  const int par, const int delta) {
  int w;
  do {
    w = intSymmWalk(r, delta);
  } while(par + w <= 0);
  return w;
}

static void changeIntPositive(const gsl_rng *r, intparameters *p, 
			      int i, int width) {
  int par=getIntParameter(p, i);
  int w=intSymmWalkPositive(r, par, width);
  setIntProposal(p, i, par+w);
}

#ifndef PROP
/* Propose parameters */
void proposeIntPositiveAB(const gsl_rng *r, intparameters *p, 
			 int i0, int steps) {
  int i, n=i0+steps;
  for(i=i0; i<n; i++) {
    changeIntPositive(r, p, i,  getIntWidth(p,i));
  }
}

/* Samples between the current values of the parameters and the maxima
   or the minima and the currrent values. */
void proposeIntIntervalAB(const gsl_rng *r, intparameters *p, 
			  int i0, int steps) {
  int i, n=i0+steps;
  for(i=i0; i<n; i++) {
    int par=getIntParameter(p, i);
    int newPar;
    int up=gsl_rng_uniform(r)<0.5;
    if(up) {
      int max=getIntMax(p, i);
      newPar=intUniform(r, par+1, max-(par+1)+1);
      setIntProposal(p, i, newPar);
      /* fprintf(OUT, ""); */
    }
    else {
      int min=getIntMin(p, i);
      newPar=intUniform(r, min, (par-min)+1);
      setIntProposal(p, i, newPar);
    }
  }
}

/* Samples between the current values of the parameters and the maxima
   or the minima and the currrent values. */
void proposeIntShiftInterval(const gsl_rng *r, intparameters *p, 
			     int i0, int nIPar) {
  /* Moves any parameter starting from i0... probably usually i0=0 */
  /* int par=getIntParameter(p, i0); */
  int newPar;
  int min=(i0>0)?getIntParameter(p,i0-1):getIntMin(p,i0);
  int max =(i0+1<nIPar)?getIntParameter(p,i0+1):getIntMax(p, i0);
  newPar=intUniform(r, min+1, max-(min+1)-1);
  /* fprintf(OUT, "moveK(): i0=%i, min=%i, max=%i, old=%i, new=%i\n", */
  /* 	  i0, min+1, max-(min+1)-1, par, newPar); */
  setIntProposal(p, i0, newPar);
  /* int up=gsl_rng_uniform(r)<0.5; */
  /* if(up) { */
  /*   /\* Sample between current value of this and next parameter or */
  /*      between current value and maximal value *\/ */
  /*   int max=(i0+1<nIPar)?getIntParameter(p,i0+1):getIntMax(p, i0); */
  /*   /\* fprintf(OUT, "k%i = %i, k%i = %i, max=%i\n", i0, par,i0+1, getIntParameter(p,i0+1),max); *\/ */
  /*   if(par+1<max-(par+1)+1) newPar=intUniform(r, par+1, max-(par+1)+1);  */
  /*   else newPar=par; */
  /*   fprintf(OUT, "k%i: up: oldK = %i, newK = %i, upper = %i\n", i0, par, newPar, max); */
  /*   setIntProposal(p, i0,newPar);  */
  /* } else {  */
  /*   /\* Sample between preceding and current value or between minimum */
  /*      and current value. *\/ */
  /*   int min=(i0-1>=0)?getIntParameter(p,i0-1):getIntMin(p, i0);  */
  /*   if(min<(par-min)) newPar=intUniform(r, min, (par-min));  */
  /*   else newPar=par; */
  /*   fprintf(OUT, "k%i, down: oldK = %i, newK = %i,  lower = %i\n", i0, par, newPar, min); */
  /*   setIntProposal(p, i0, newPar);  */
  /* }  */
}


static void changeIntConstrained(const gsl_rng *r, intparameters *p, 
				 int i) {
  int par=getIntParameter(p, i), w;
  int width = getIntWidth(p,i), min=getIntMin(p,i), max=getIntMax(p,i);

  do {
    w = intSymmWalk(r, width);
/*     w = intUpDown(r, width); */
    /* fprintf(OUT, "changeIntConstrained(): width = %i, par =%i, w = %i\n",  */
    /* 	    width, par, w); */
  } while( (par + w < min) || (par + w > max) );
  setIntProposal(p, i, par+w);

/*   do { */
/*     sample = intUniform(r, par-width,2*width+1); */
/*   } while( (sample < min) || (sample > max) ); */
/*   setIntProposal(p, i, sample); */
}

/* Propose parameters */
void proposeIntConstrainedAB(const gsl_rng *r, intparameters *p, int i0, int steps) {
  int i;
  for(i=i0; i<i0+steps; i++) {
    changeIntConstrained(r, p, i);
  }
}


/* check if proposal is positive */
int isIntProposalPositive(const intparameters *p, int nPar) {
  int i;

  for(i=0; i<nPar; i++) {
    int prop = getIntProposal(p, i);
    if( prop < EPS) {
/*       fprintf(ERR, "negative proposal found:\n"); */
      return 0;
    }
  }
  return 1;
}
#endif
/******************* End: int parameters ********************/
