/*
 * nw_data.c
 *
 * Read data from text files
 *
 * Ivo Siekmann, 01/07/2014
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nw_data.h"

#define ERR stderr
#define OUT stdout

/* Nach Schreiners C-Vorlesung:
 *
 *	split(0, s)	legt einen Vektor an, der auf Worte in s[] zeigt
 *	split(v, s)	verlaengert den Vektor
 *
 *	der Vektor hat minimale Laenge und NULL am Schluss
 *	s[] wird mit '\0' unterteilt
 *
 *	ist keinerlei Platz vorhanden, so ist das Resultat NULL
 *	andernfalls ist s[] moeglicherweise nicht ganz zerlegt
 */

#define	VIELE	1000
/* #define	DELIM	" \f\n\r\t\v"	/\* Trenner *\/ */
#define	DELIM	"\n\t"	/* Trenner */

#define BIGBUFSIZE 1000

#define	new(n)	    (char **) malloc((n) * sizeof(char *))
#define	renew(p, n) (char **) realloc((p), (n) * sizeof(char *))

static char ** split (char ** list, char * s){	
  char * cp, ** lp;
  unsigned lim;

  if (list) {
    /* lim bis zum Ende von list bewegen, lim ist die Anzahl der
       Eintraege */
    for (lp = list; * lp; ++ lp);
    lim = lp - list + 1;
  }
  else if ((list = new(VIELE)))
    lp = list, lim = VIELE;
  /* kein Speicher vorhanden */
  else return NULL;

  /* cp enthaelt in jedem Schritt einen Zeiger auf das naechste
     Token... und lp enthaelt jeweils die cp */
  for (cp = strtok(s, DELIM); cp; cp = strtok(NULL, DELIM)) {
    if (lp - list == lim - 1) {
      if (! (lp = renew(list, lim + VIELE))) {
	list[lim-1] = NULL;
	return list;
      }
      list = lp, lp = list + lim - 1, lim += VIELE;
    }
    * lp ++ = cp;
  }
  /* Letzter Eintrag von lp ist NULL */
  * lp ++ = NULL;
  return renew(list, (lp - list));
}

static char * strnsave (const char * s, size_t n) {
  char * p = malloc((strlen(s) + n + 1) * sizeof(char));

  if (! p) {
    fprintf(stderr, "strsave: no room");
    exit(1);
  }

  return strcpy(p, s);
}

static char * strsave (const char * s)
{
	return strnsave(s, 0);
}


/* Converts string s to tokenizes line nw_tokline. Return NULL if this failed. */
nw_tokline * nw_data_strToTokline(const char *s) {
  /* String MUST be copied - otherwise only the last row is actually saved! */
  char *newS=strsave(s);
  char **tokenizedline=split(NULL, newS);
  nw_tokline *tokline;
  
  /* fprintf(stderr, "nw_data_strToTokline(): Tokenizing string %s\n", s); */
  if(!tokenizedline) return NULL;

  tokline=malloc(sizeof(*tokline));
  if(!tokline) return NULL;

  tokline -> s=tokenizedline;

  /* nw_data_printTokline(stdout, tokline); */

  /* fprintf(stderr, "\n"); */

  return tokline;
}

/* Reads file fp line by line and save tokenized lines in tl. Returns number of lines successfully read. */
int nw_data_fpToTokline(FILE *fp, nw_tokline **tl) {
  /*Text-Datei wird zunaechst in "Tokens" zerlegt*/
  char buf[BIGBUFSIZE];
  /*Anzahl der Zeilen*/
  unsigned lno = 0;

  /* setvbuf(fp, buf, _IOLBF, BIGBUFSIZE); */
  
  while (fgets(buf, sizeof buf, fp)) {
    nw_tokline *l;
    if (buf[strlen(buf)-1] != '\n') {
      fprintf(stderr, "line %u too long:\n", lno);
      fprintf(stderr, "%s\n", buf);
      exit(1);
    }
    /* fprintf(stderr, "nw_data_fpToTokline(): Tokenizing line %s\n", buf ); */
    l=nw_data_strToTokline(buf);
    if(!l) {
      fprintf(stderr, "nw_data_fpToTokline(): Failed reading from file\n");
    }
    tl[lno++]=l;
  }
  return lno;
}

/* Count the number of entries in tl */
int nw_data_nCol(const nw_tokline *tl) {
  char **cp;
  int nCol=0;
  for(cp=tl->s, nCol=0; *cp; ++cp, nCol++);

  return nCol;
}

/* Print tokenized line */
void nw_data_printTokline(FILE *fp, const nw_tokline *tl) {
  char **cp;

  for(cp=tl->s; *cp; ++cp) {
    fprintf(fp, "%s\t", *cp);
  }
  fprintf(fp, "\n");
}


char *nw_data_getRowJ(const nw_tokline *tl, int j) {
  int k;
  char **cp;
  for(cp=tl->s, k=0; *cp && k<j; ++cp, k++) {
    fprintf(stderr, "%s\t", *cp);
  }
  fprintf(stderr, "\n");
  return *cp;
}

/* Return string at row i, column j */
char *nw_data_getStrIJ(const nw_tokline **tm, int i, int j) {
  const nw_tokline *tl;
  char **cp;
  
  tl=tm[i];
  /* fprintf(ERR,"nw_data_getStrIJ(): tm[%i]=%p\n", i, tl); */
  cp=tl->s;
  /* fprintf(ERR,"nw_data_getStrIJ(): s(%i,%i)=%p(%s)\n", i, j, cp, cp[j]); */
  return cp[j];
}

/* Print tokenized lines */
void nw_data_printTokmatrix(FILE *fp, const nw_tokline **tm, int n) {
  int i;

  for(i=0; i<n; i++) {
    const nw_tokline *tl=tm[i];
    /* fprintf(stderr, "nw_data_printTokmatrix(): Printing line %i\n", i); */
    nw_data_printTokline(fp, tl);
  }
}

/* If string s is in nw_tokline*/
int nw_data_findString(const nw_tokline *tl, const char *s) {
  char **cp;
  int k;
  for(cp=tl->s, k=0; *cp; ++cp, k++) {
    /* If string s found return index k */
    /* fprintf(ERR, "nw_data_findString(): %i: %s (%p)\n", k, *cp, *cp); */
    if(!strcmp(*cp, s)) return k;
  }
  
  return -1;
}

/* If string s is 'row header' of one of the rows of tm return index. Otherwise return -1 */
int nw_data_findRowHeader(const nw_tokline **tm, const char *s, int N) {
  int i;

  /* fprintf(ERR, "nw_data_findRowHeader(): Finding row header %s...\n", s); */
  for(i=0; i<N; i++) {
    char *strIJ=nw_data_getStrIJ(tm, i,0);
    /* fprintf(ERR, "nw_data_findRowHeader(): %i - %s\n", i, strIJ); */
    /* Ignore empty lines */
    if(strIJ){
      if(!strcmp(strIJ, s)) return i;
    }
  }

  return -1;
}

/* If string s is 'column header' of one of the rows of tm return index. Otherwise return -1 */
int nw_data_findColumnHeader(const nw_tokline **tm, const char *s) {
  /* fprintf(ERR, "nw_data_findRowHeader(): Finding column header %s...\n", s); */
  return nw_data_findString(tm[0], s);
}

/* Select by row header rh and column header ch. Return NULL if either of them is not found. N is maximum number of lines. */
char *nw_data_selectRowCol(const nw_tokline **tm, 
			   const char *rh, const char *ch, int N) {
  int i=nw_data_findRowHeader(tm, rh, N), j=nw_data_findColumnHeader(tm, ch);

  if( (i<0) || (j<0)) return NULL;

  /* fprintf(ERR, "nw_data_selectRowCol(): Selecting string at (%s, %s)...\n", */
  /* 	  rh, ch); */

  return nw_data_getStrIJ(tm, i, j);
}

/* Select by row index i and column header ch. Return NULL if either
   of them is not found. N is maximum number of lines. */
char *nw_data_selectRowIndCol(const nw_tokline **tm, 
			      int i, const char *ch, int N) {
  int j=nw_data_findColumnHeader(tm, ch);

  if( (i<0) || (j<0)) return NULL;

  /* fprintf(ERR, "nw_data_selectRowCol(): Selecting string at (%d, %s)...\n", */
  /* 	  i, ch); */

  return nw_data_getStrIJ(tm, i, j);
}

/* Select string by row and column header and convert it to double */
double nw_data_selectDoubleByRowCol(const nw_tokline **tm, 
				    const char *rh, const char *ch, int N) {
  char *s=nw_data_selectRowCol(tm, rh, ch,  N);

  double d;
  if(!s) fprintf(ERR, 
		 "nw_data_selectDoubleByRowCol(): Failed  to select (%s,%s)\n",
		 rh,ch), exit(1);

  d  =strtod(s, NULL);
  if(!d) fprintf(ERR, "nw_data_selectDoubleByRowCol(): Failed to convert %s to double\n", s), exit(1);

  return d;
}

/* Select string by row and column header and convert it to double */
double nw_data_selectDoubleByRowIndCol(const nw_tokline **tm, 
				       int i, const char *ch, int N) {
  char *s=nw_data_selectRowIndCol(tm, i, ch,  N);

  double d;
  if(!s) fprintf(ERR, 
		 "nw_data_selectDoubleByRowCol(): Failed  to select (%d,%s)\n",
		 i,ch), exit(1);

  d  =strtod(s, NULL);
  if(!d) fprintf(ERR, "nw_data_selectDoubleByRowCol(): Failed to convert %s to double\n", s), exit(1);

  return d;
}

/* Select multiple doubles by row and multiple nCols column header
   from a matrix with N rows.*/
double *nw_data_selectDoublesByRowColumns(const nw_tokline **tm, 
					  const char *rh, const char *ch[], 
					  int nCols, int N) {
  double *v=malloc(nCols*sizeof(*v));
  int i;
  
  for(i=0; i<nCols; i++) {
    double d=nw_data_selectDoubleByRowCol(tm, rh, ch[i], N);
    v[i]=d;
  }
  return v;
}

/* Return N strings in i-th column, starting from row j0 */
char **nw_data_getColumn(const nw_tokline **tm, int i, int j0,  int N) {
  char **out=malloc(N*sizeof(*out));
  int j, k;

  /* Does string need to be COPIED? Or is it ok to just return a
     pointer? */
  for(j=j0, k=0; j<j0+N; j++, k++) {
    char *s=nw_data_getStrIJ(tm, j, i);
    out[k]=s;
  }
  
  return out;
}


/* Return N strings from first column where s is the 'header' */
char **nw_data_getColumnByHeader(const nw_tokline **tm, const char *s, int N) {
  int i=nw_data_findColumnHeader(tm, s);

  if(i<0) return NULL;

  return nw_data_getColumn(tm, i, 1,  N);
}

/* SHOULD return index until which conversion was successful */
/* Convert N strings in i-th column, starting from row j0 to doubles */
int nw_data_columnToDouble(const nw_tokline **tm, int i, int j0,  int N,
			    double *v) {
  char **col=nw_data_getColumn(tm, i, j0,  N);
  int j;

  for(j=0; j<N; j++) {
    double d;
    /* fprintf(ERR, "nw_data_columnToDouble(): j=%i: %s\n", j, col[j]); */
    d=strtod(col[j], NULL);
    if(!d) { return j;}
    v[j]=d;
  }

  return N;
}

/* Converts columns nCol columns starting from index col0 to
   double. In each column the first entry is j0.
*/
int nw_data_convertColumnsToDouble(const nw_tokline **tm, int col0, int nCol,
				   int j0,  int N, double *v) {
  int i;
  int nDoubles=0;
  for(i=col0, nDoubles=0; i<col0+nCol; i++) {
    int nTrace;
    nTrace=nw_data_columnToDouble(tm,
				  i, j0,  N, v+nDoubles);
    nDoubles+=nTrace;
  }
  return nDoubles;
}

/* SHOULD return index until which conversion was successful */
/* Return N strings from first column where s is the 'header' */
int nw_data_columnToDoubleByHeader(const nw_tokline **tm, const char *s, int N, double *v) {
  int i=nw_data_findColumnHeader(tm, s);

  if(i<0) return -1;

  return nw_data_columnToDouble(tm, i, 1, N, v);
}

/* Return N strings from i-th row, starting from row j0 */
char **nw_data_getRow(const nw_tokline **tm, int i, int j0,  int N) {
  char **out=malloc(N*sizeof(*out));
  int j, k;

  for(j=j0, k=0; j<j0+N; j++, k++) {
    char *s=nw_data_getStrIJ(tm, i, j);
    out[k]=s;
  }
  
  return out;
}
