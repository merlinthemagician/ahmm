#ifndef SAMPLE
#define SAMPLE
/*
 * mp_sample.h
 *
 * MCMC parameter data type
 * 
 * 
 *
 * Ivo Siekmann, 28/06/2013
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */

#include <stdio.h>
#include <stdlib.h> 

/* #include <gsl/gsl_rng.h> */
/* #include <gsl/gsl_matrix.h> */

/* typedef double samp_t; */
typedef int samp_t;

/* Data structure for samples of generic type */
typedef struct mp_sample{
  size_t n;
  samp_t *data;
} mp_sample;

typedef mp_sample *mp_sample_set[];

/* allocates memory for n parameters */
mp_sample *mp_sample_alloc(int n);

/* Read samples from file and save to mp_sample_set */
void mp_sample_fscanf(FILE *fp, char *fmt, mp_sample_set s);

double *mp_sample_getDoubles(const mp_sample* sample);
#endif
