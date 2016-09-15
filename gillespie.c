/*
 * gillespie.c
 *
 * Simulation of an N-state continuous Markov chain model with the
 * Gillespie algorithm (Gillespie, Journal of Physical Chemistry,
 * Vol. 81, 1977).
 *
 * Ivo Siekmann, April 17th 2009
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */

#include <stdio.h>

#define OUT stdout
#define ERR stderr


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>


/* active state */
static int state = 3, newState=3;

/* current time*/
static double t = 0;

/* sampling time */
static double tSamp=0;


/* waiting time */
static double waitT = 0;

/* sum of the concentrations times rates */
static double a0 = 0;

/* random number generator */
static gsl_rng * r;

/* File name of model file */
static char* modelFn;

/* Number of total model states and open states. */
static int nStates, nOpen;

/* Infinitesimal generator */
static double **Qrates;

/* End time */
static double tend;

/* Initial value of random number generator */
static int seed=42;

/* Standard deviation of normally distributed noise */
static double noiseIntensity;

/* Sampling interval */
static double tau=0.05;

/* Open conductance */
static double openConductance=-20;


/* Compute next state of the Markov chain */
void next() {
  double rT;
  int j;

  double rThresh = 0, sumOfRates = 0;

  /* Generating uniformly distributed random numbers for computation
   * of: 1) waiting time and 2) transition */
  rT = gsl_rng_uniform (r);

  a0= -Qrates[state][state];

  /*Compute waiting time */
  waitT = log(1.0/rT)/(a0);

  rThresh = a0*gsl_rng_uniform(r);

  /* Calculate which state transition has occured during the waiting
     time, skip index of state*/ 
  for(j=0; (j<nStates) && (rThresh > sumOfRates); j++) {
    /* Ignore Q_ii */
    if(j == state) j++;
    sumOfRates += Qrates[state][j]; 
  }

  /* current state -> j */
  newState = j-1;
}

static double noiseIntensity;

static int verbose=0;

/*Prints data - possibly after processing*/
void printData(FILE *f, double t, int state) {
  double curr = (state>=(nStates-nOpen))?openConductance:0, noise;

  noise =gsl_ran_gaussian(r, noiseIntensity);

  /* Print either:
     time, original record, noisy record, state 
     or:
     time, noisy record
  */
  if (verbose) 
    fprintf(f,"%f\t%f\t%f\t%i\n", t, curr, curr+noise,state);
  else
    fprintf(f,"%f\t%f\n", t, curr+noise);
}

/* Output data to file f */
void foutput(FILE *f) {
  printData(f, t, state);
}

/* Output data */
void out() {
  foutput(OUT);
}

/* Output data in the "sampling interval" tau */
void fsampleout(FILE *f, double tau) {
  for(;tSamp<t+waitT; tSamp+=tau) {
    printData(f, tSamp, state);
  }
}

/* Output data in the "sampling interval" tau */
void sampleout(double tau) {
  fsampleout(OUT, tau);
}

/* Run till tend */
void frun(FILE *f, double tend) {
  do{
    foutput(f);
    next();
    t += waitT;
    state = newState;
  } while(t/* +waitT */<tend);
}

/* Run till tend */
void run(double tend) {
  frun(OUT, tend);
}

/* Run giving sampled output with sampling interval tau till tend */
void fsamprun(FILE *f, double tau, double tend) {
  do {
    t += waitT;
    state = newState;
    next();
    fsampleout(f,tau);
  } while(tSamp<tend);
}

/* Run giving sampled output with sampling interval tau until tend */
void samprun(double tau, double tend) {
  fsamprun(OUT,tau, tend);
}

void makeConservative() {
  int i, j;
  double sumOfRates = 0;
  for(i=0; i<nStates; i++) {
    sumOfRates = 0;
    for(j=0; j<i; j++) {
      sumOfRates += Qrates[i][j];
    }
    for(j=i+1; j<nStates; j++) {
      sumOfRates += Qrates[i][j];
    }
    Qrates[i][i] = -sumOfRates;
  }
}

void processArgs(int argc, char **argv){
  int i,j;
  int nArgs=5, nArgMax=9;
  const char *usage="modelfile nStates nOpen tend [seed] [noiseIntensity] [samplingInterval] [openConductance]";
  char * endptr;
  FILE *fp;
  
  if( (argc<nArgs) || (argc>nArgMax)) {
    fprintf(ERR, "Usage:\n%s %s\n", argv[0], usage);
    exit(1);
  }

  modelFn=argv[1];

  nStates=strtol(argv[2], &endptr, 10);
  if(*endptr != '\0') {
    fprintf( ERR, "Error converting %s to integer\n", argv[2]);
    exit(1);
  }

  double (*qr)[nStates]=malloc(sizeof (double[nStates][nStates]));
  fp=fopen(modelFn,"r");
  
  io_readMatrix(fp, nStates, nStates, qr);
  /* copy to Qrates */
  Qrates=malloc(nStates*sizeof(double*));
  for(i=0; i<nStates; i++) {
    Qrates[i] = malloc(nStates*sizeof(double));
    for(j=0; j<nStates; j++) {
      Qrates[i][j]=qr[i][j];
    }
  }
  fclose(fp);

  nOpen=strtol(argv[3], &endptr, 10);
  if(*endptr != '\0') {
    fprintf( ERR, "Error converting %s to integer\n", argv[3]);
    exit(1);
  }

  tend=strtod(argv[4], &endptr);
  if(*endptr != '\0') {
    fprintf( ERR, "Error converting %s to double\n", argv[4]);
    exit(1);
  }

  if(argc>nArgs) {
    seed=strtol(argv[5], &endptr, 10);
    if(*endptr != '\0') {
      fprintf( ERR, "Error converting %s to integer\n", argv[5]);
      exit(1);
    }
  }

  if(argc>nArgs+1) {
    noiseIntensity=strtod(argv[6], &endptr);
    if(*endptr != '\0') {
      fprintf( ERR, "Error converting %s to double\n", argv[6]);
      exit(1);
    }
  }

  if(argc>nArgs+2) {
    tau=strtod(argv[7], &endptr);
    if(*endptr != '\0') {
      fprintf( ERR, "Error converting %s to double\n", argv[7]);
      exit(1);
    }
  }

  if(argc>nArgs+3) {
    openConductance=strtod(argv[8], &endptr);
    if(*endptr != '\0') {
      fprintf( ERR, "Error converting %s to double\n", argv[8]);
      exit(1);
    }
  }
}


int main(int argc, char **argv) {
  double tau =0.05, tend = 2000;
  long int seed = 42;
  int i,j;

  noiseIntensity = sqrt(7.5);
  processArgs(argc, argv);
    
  r = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(r, seed);
  state=gsl_rng_uniform_int(r,nStates);
  newState=state;

  makeConservative();

  fprintf(ERR, "\nMatrix des Markov-Modells:\n");
  for(i=0; i<nStates; i++) {
    for(j=0; j<nStates; j++) {
      fprintf(ERR, "%g\t", Qrates[i][j]);
    }
    fprintf(ERR,"\n");
  }
  fprintf(ERR,"\n");
  
  fprintf(ERR, "Starting in state %i\n", state);

  samprun(tau,tend);

  gsl_rng_free(r);

  return 0;
}
