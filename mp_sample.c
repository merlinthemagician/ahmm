/*
 * mp_sample.c
 *
 * Data type for samples
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

#include "mp_sample.h"

/* allocates memory for n parameters */
mp_sample* mp_sample_alloc(int n) {
  mp_sample* s=malloc(sizeof(mp_sample));
  s->n=n;
  s->data=malloc(n*sizeof(samp_t));

  return s;
}

/* Read samples from file and save to mp_sample_set */
/* void mp_sample_fscanf(FILE *fp, char *fmt, mp_sample_set s); */

/* double *mp_sample_getDoubles(const mp_sample* sample); */
