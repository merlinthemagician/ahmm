/*
 *  matrixIO.c
 *
 *
 * Ausgaberoutinen fuer Matrizen
 *
 * Ivo Siekmann, 31.10.2006
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */
#include "matrixIO.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>

/* for LINE_MAX */
#include <limits.h>

#define OUT stdout
#define ERR stderr

/*Ausgabe eines gsl_vector*/
void vprintfmt(FILE * stream, gsl_vector * v, const char* fmt) {
  int i;
  const char del='\t';
  for(i=0; i<v->size-1; i++) {
    fprintf(stream,fmt,gsl_vector_get(v,i));
    fputc(del, stream);
  }
  fprintf(stream,"%f\n",gsl_vector_get(v,i));
}

/*Ausgabe eines gsl_vector*/
void vprint(FILE * stream, gsl_vector * v) {
  vprintfmt(stream, v, "%f");
}


/*
 * Ausgabe eines Zeitschritts nach f
 */
void
output(FILE* f, double t, gsl_vector* y)
{
  fprintf(f,"%f ",t);
   vprint(f,y);
}


/*Ausgabe eines Zeitschritts*/
void printStep(FILE * stream, double t, double hx, gsl_vector * u) {
  int i;

  for(i=0; i<u->size; i++) {
    fprintf(stream, "%f\t", t);
    fprintf(stream, "%f\t", i*hx);
    fprintf(stream, "%f\n", gsl_vector_get(u,i));
  }
  fprintf(stream, "\n");
}

/* Ausgabe einer gsl_matrix in das File f. Zeilentrenner rd und
   Spaltentrenner cd sowie Formatstring fmt */
void printMatf(FILE *f, const gsl_matrix *m,
	       char *rd,char *cd, char *fmt) {
  int i, j;
  for(i=0; i<m->size1; i++) {
    for(j=0; j<m->size2-1; j++) {
      fprintf(f, fmt, gsl_matrix_get(m, i,j));
      fprintf(f, "%s", cd);
    }
    fprintf(f, fmt, gsl_matrix_get(m, i,j));
    fprintf(f, "%s", rd);
  }
}

/*Ausgabe einer gsl_matrix in das File f*/
void printMat(FILE *f, const gsl_matrix *m) {
  printMatf(f, m, "\n", "\t", "%f");
}

/*Umwandlung von Ziffern in entsprechende ASCII-Chars*/
char i2c(int i) {

  switch(i) {
  case 0: {
    fprintf(stderr,"Zu wenig Stellen fuer Nummerierung des Dateinamen\n");
    exit(1);
   }
  case 1: return '1';
  case 2: return '2';
  case 3: return '3';
  case 4: return '4';
  case 5: return '5';
  case 6: return '6';
  case 7: return '7';
  case 8: return '8';
  case 9: return '9';
  default: {
    fprintf(stderr, "Nur 9 Stellen fuer Nummerierung unterstuetzt\n");
    exit(1);
  }
  }
}

/*Konstruiere einen Dateinamen zur Ausgabe von nummerierten Dateien:
 Praefix pre, Suffix suf, Nummer n, verfuegbare Ziffern digits*/ 
char* numberedFN(const char* pre, const char* suf, int n, int digits) {

  int len = strlen(pre) + digits + strlen(".dat")+1;
  char* outname = (char*)malloc(len*sizeof(char));
  char mid[] = "%00i";
  char* fmt;

  /*Setze in mid die Laenge ein*/
  mid[2] = i2c(digits);

  /*Konstruiere Formatierungsstring*/
  len = strlen(pre)+strlen(mid) + strlen(suf)+1;
  fmt = (char*)malloc(len*sizeof(char));
  strcpy(fmt, pre);

  fmt = strcat(fmt, mid);
  fmt = strcat(fmt, suf);

  /*Konstruiere Dateinamen*/
  sprintf(outname, fmt , n);

  free(mid);
  free(fmt);

  return outname;
}

/* Ausgabe des Zeitschritts i
 * in Datei mit Namen "pre"+Nummerierung
 * digits darf nicht groesser als 9 sein!
 */ 

void matToFile(const gsl_matrix *mat, const char* pre, int i, int digits) {

  char* suf = ".dat";
  char* outname;/*  = (char*)malloc(len*sizeof(char)); */
  FILE *f;

/*   printf("In matToFile():\n"); */

  /*Konstruiere Dateinamen*/
  outname = numberedFN(pre, suf, i, digits);
/*   printf(outname); */
  /*Ausgabe in Datei*/
  f = fopen(outname, "w");
  printMat(f, mat);
  fflush(f);
  fclose(f);
}

/************************************************************************/
/*2D-Ausgabe*/
/************************************************************************/

/*Erzeugt einen Dateinamen mit Praefix pre*/
/* char * genfilename(char* pre, char* suf, int eqn, double t) { */
/*   char * s = malloc(strlen(pre)+strlen("u0t")+13+strlen(suf)+1); */
/*   char *fmt = malloc(strlen(pre)+strlen("u%it%09.3f")+1); */

/*   strcpy(fmt,pre); */
/*   strcat(fmt,"u%it%09.3f"); */
/*   sprintf(s, fmt, eqn, t); */
/*   strcat(s,suf); */

/* /\*   free(fmt); *\/ */

/*   return s; */
/* } */

/*Schreibt den Kopf fuer die Ausgabedatei*/
void writeHead(FILE* f, char *s) {
  
}

/*
 * Schreibt Farbwert fuer eine binaere PNM-Datei: Der double-Wert wird
 * mit s skaliert.
 */
void writeColourBin(FILE *f, double d, double s) {
  int colour;

  colour = (int)(d*s);
  if(colour > 255) colour = 255;
  putc(colour,f);
}

/*
 * Schreibt Farbwert fuer eine textbasierte PNM-Datei: Der double-Wert
 * wird mit s skaliert.
 */
void writeColourText(FILE *f, double d, double s) {
  int colour;

  colour = (int)(d*s);
  if(colour > 255) colour = 255;
  fprintf(f, "%i ",colour);
}

/*
 * Schreibt Farbwert fuer eine binaere PNM-Datei: Der double-Wert wird
 * mit s skaliert und danach invertiert.
 */
void writeColourBinInv(FILE *f, double d, double s) {
  int colour;

  colour = (int)(d*s);
  if(colour > 255) colour = 255;
  putc(255-colour,f);
}

/*
 * Schreibt Farbwert fuer eine textbasierte PNM-Datei: Der double-Wert
 * wird mit s skaliert und anschliessend invertiert.
 */
void writeColourTextInv(FILE *f, double d, double s) {
  int colour;

  colour = (int)(d*s);
  if(colour > 255) colour = 255;
  fprintf(f, "%i ", 255-colour);
}


/*
 * Gibt Loesung von drei Spezies als RGB-Datei aus: fn ist ein
 * Dateiname mit Suffix '.ppm'. Die Daten sind im gsl_vector u
 * gespeichert. Es werden die Indizes indRot, indGruen und indBlau
 * verwendet, um den Wert fuer die "rote, die "gruene" und die "blaue"
 * Population zu lesen. Die Werte werden skaliert mit den Faktoren
 * sRot, sGruen und sBlau. Die Raumdimensionen des zweidimensionalen
 * Gitters sind dimx und dimy.
 */
void matToPPM(char *fn, gsl_vector_view *uView, 
	      int indRot, int indGruen, int indBlau,
	      double sRot, double sGruen, double sBlau, 
	      int dimx, int dimy) {
  FILE *f;
  double d;
  int i;

  f = fopen(fn,"w");

  /*Schreibe Dateikopf*/
  fprintf(f, "P6\n");

  /*Gitterdimensionen*/
  fprintf(f, "%i %i\n",dimx, dimy);
  fprintf(f,"255\n");

  /*Schreiben der Farbdaten*/
  for(i=0; i<dimx;i++) {
    d = gsl_vector_get(&uView[i].vector, indRot);
    writeColourBin(f,d,sRot);
    d = gsl_vector_get(&uView[i].vector, indGruen);
    writeColourBin(f,d,sGruen);
    d = gsl_vector_get(&uView[i].vector, indBlau);
    writeColourBin(f,d,sBlau);
  }
  fflush(f);
  fclose(f);
}


/*
 * Gibt Loesung von drei Spezies als RGB-Datei aus: fn ist ein
 * Dateiname mit Suffix '.ppm'. Die Daten sind im gsl_vector u
 * gespeichert. Es werden die Indizes indRot, indGruen und indBlau
 * verwendet, um den Wert fuer die "rote, die "gruene" und die "blaue"
 * Population zu lesen. Die Werte werden skaliert mit den Faktoren
 * sRot, sGruen und sBlau. Die Raumdimensionen des zweidimensionalen
 * Gitters sind dimx und dimy.
 */
void printMatPPM(char *fn, gsl_matrix_view *uEqnView, 
	      int indRot, int indGruen, int indBlau,
	      double sRot, double sGruen, double sBlau) {
  FILE *f;
  double d;
  int i,j;

  f = fopen(fn,"w");

  /*Schreibe Dateikopf*/
  fprintf(f, "P6\n");

  /*Gitterdimensionen*/
  /* Use C99 flag for output */
  fprintf(f, "%zu %zu\n",uEqnView[0].matrix.size1, uEqnView[0].matrix.size2);
  fprintf(f,"255\n");

  /*Schreiben der Farbdaten*/
  for(i=0; i<uEqnView[0].matrix.size1;i++) {
    for(j=0; j<uEqnView[0].matrix.size2;j++) {
      d = gsl_matrix_get(&uEqnView[indRot].matrix, j, i);
      writeColourBin(f,d,sRot);
      d = gsl_matrix_get(&uEqnView[indGruen].matrix, j, i);
      writeColourBin(f,d,sGruen);
      d = gsl_matrix_get(&uEqnView[indBlau].matrix, j, i);
      writeColourBin(f,d,sBlau);
    }
  }
  fflush(f);
  fclose(f);
}


void writePGMHead(FILE *f, int dimx, int dimy) {
#ifdef BINPNM
  fprintf(f, "P5\n");
#else
  fprintf(f, "P2\n");
#endif
  /*Gitterdimensionen*/
  fprintf(f, "%i %i\n",dimx, dimy);

  fprintf(f,"255\n");
}


/*Gibt Loesung als Grauwerte in PGM-Format aus*/
void matToPGM(char *fn, gsl_vector_view *uView,
	      int ind, double s, int dimx, int dimy) { 

  FILE *f;
  int i;
  double d;

  f = fopen(fn,"w");

  writePGMHead(f, dimx, dimy);
  /*Farbe wird invertiert*/
  for(i=0; i<dimx*dimy;i++) {
      d = gsl_vector_get(&uView[i].vector, ind);
#ifdef BINPNM
      writeColourBinInv(f,d,s);
#else
      writeColourTextInv(f,d,s);
#endif
#ifndef BINPNM
      if(i%dimy==0)
	fprintf(f,"\n");
#endif
  }
  fflush(f);
  fclose(f);
}


/*Ausgabe einer gsl_matrix in das File f*/
void printMatPGM(FILE *f, const gsl_matrix *m, double s) {
  int i, j;
  double d;

  /*X/Y vertauscht !*/
  /* Vertauschen scheint noetig zu sein*/
  writePGMHead(f, m->size2, m->size1);
  for(i=0; i<m->size1; i++) {
    for(j=0; j<m->size2; j++) {
      /*X/Y vertauscht !*/
      /* Vertauschen scheint noetig zu sein*/
      d = gsl_matrix_get(m, j, i);
#ifdef BINPNM
      writeColourBinInv(f,d,s);
#else
      writeColourTextInv(f,d,s);
#endif
    }
#ifndef BINPNM
    fprintf(f,"\n");
#endif
  }
}


/*Gibt Loesung als double-Werte aus*/
void matToDat(char *fn, gsl_vector_view *uView,
	      int ind, int dimx, int dimy) { 

  FILE *f;
  int i;
  double d;

  f = fopen(fn,"w");

  for(i=0; i<dimx*dimy;i++) {
      d = gsl_vector_get(&uView[i].vector, ind);
      fprintf(f, "%f\t", d);
      if(i%dimy==0)
	fprintf(f,"\n");
  }
  fflush(f);
  fclose(f);
}

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

#define	VIELE	4000
#define	DELIM	" \f\n\r\t\v"	/* Trenner */

#define	new(n)	    (char **) malloc((n) * sizeof(char *))
#define	renew(p, n) (char **) realloc((p), (n) * sizeof(char *))

char ** split (char ** list, char * s){	
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

char * strnsave (const char * s, size_t n) {
  char * p = malloc((strlen(s) + n + 1) * sizeof(char));

  if (! p) {
    fprintf(stderr, "strsave: no room");
    exit(1);
  }

  return strcpy(p, s);
}

char * strsave (const char * s)
{
	return strnsave(s, 0);
}

double * dblnsave (const double * s, size_t n) {
  double * p = malloc( n * sizeof(double));
  int i;

  if (! p) {
    fprintf(stderr, "dblsave: no room");
    exit(1);
  }
  for(i=0; i<n; i++) {
    p[i]=s[i];
  }
  return p;
}

double * dblsave (const double * s){
  return dblnsave(s, 0);
}

/* Counts fields in fp delimited by DELIM */
int countFields(FILE *fp) {
  /*count columns and delete first line */
  char line[LINE_MAX], *cp;
  int count=0;
  if(fgets(line, sizeof line, fp) == NULL) {
    fprintf(ERR, "countFields(): Error reading file\n");
  }
  for (cp = strtok(line, DELIM); cp ; cp = strtok(NULL, DELIM)) count++;
 
  rewind(fp);
  return count; 
}

/*Reads a dynamically allocated double matrix. Returns number of
  rows. */
int io_readMatrix(FILE *fp, int nRows, int nCols, double (*m)[nCols]) {
  int k=0, l=0;
  double d;
  for(k=0; k<nRows; k++) {
    for(l=0; l<nCols; l++) {
      int err=fscanf(fp,"%lf", &d);
      if( err==EOF) {
	fprintf(ERR, "readMatrix(): File ended in line %i, before M=%i\n",
		k,  nRows);
	return k;
      }
      if(!err) {
	fprintf(ERR, "readMatrix(): Error reading line k=%i\n", k), exit(1);
      }
      m[k][l]=d;
    }
  }
  return k;
}

#define BIGBUFSIZE 10000

/*Liest eine Matrix mit double-Werten aus einer Datei ein. Es wird
  zunaechst von einer quadratischen Matrix ausgegangen*/
gsl_matrix *readMat(FILE *fp) {
  /*Auszugebende Matrix*/
  gsl_matrix *m;
  /*Anzahl der enthaltenen Zahlen und Zeilen*/
  int i=0;

  /*Text-Datei wird zunaechst in "Tokens" zerlegt*/
  char ** lp, ** list = split(NULL, "");
  char buf[BIGBUFSIZE];
  /*Anzahl der Zeilen*/
  unsigned lno = 0;

  setvbuf(fp, NULL/*buf*/, _IOLBF, BIGBUFSIZE);

  if (! list)
    fputs("no room\n", stderr), exit(1);

  while (fgets(buf, sizeof buf, fp))
    {	++ lno;
      if (buf[strlen(buf)-1] != '\n') {
	fprintf(stderr, "line %u too long\n", lno);
	exit(1);
      }
      list = split(list, strsave(buf));
    }

  fprintf(stdout, "Zeilen: %i\n", lno);
  /*Verarbeite Strings zu doubles*/
  m = gsl_matrix_alloc(lno, lno);
  for (lp = list, i=0; * lp; ++ lp, i++) {
    gsl_matrix_set(m, i/lno, i%lno, atof(*lp));
  }

  return m;
}

/*Liest eine Matrix mit double-Werten aus einer Datei ein. Es wird
  zunaechst von einer quadratischen Matrix ausgegangen, Speicher muss
  reserviert sein.*/
void doubleMat(FILE *fp, double *m[], int lno) {
  /*Anzahl der enthaltenen Zahlen und Zeilen*/
  int i=0;

  /*Text-Datei wird zunaechst in "Tokens" zerlegt*/
  char ** lp, ** list = split(NULL, "");
  char buf[BIGBUFSIZE];

  setvbuf(fp, NULL/*buf*/, _IOLBF, BIGBUFSIZE);

  if (! list)
    fputs("no room\n", stderr), exit(1);

  while (fgets(buf, sizeof buf, fp))
    {/* 	++ lno; */
      if (buf[strlen(buf)-1] != '\n') {
	fprintf(stderr, "line %u too long\n", lno);
	exit(1);
      }
      list = split(list, strsave(buf));
    }

  /*Verarbeite Strings zu doubles*/
  for (lp = list, i=0; * lp; ++ lp, i++) {
    m[i/lno][i%lno] = atof(*lp);
  }
}

/*Liest eine Matrix mit double-Werten aus einer Datei ein. Es wird
  zunaechst von einer quadratischen Matrix ausgegangen, Speicher muss
  reserviert sein. KEINE Warnung bei Ueberschreiten der maximalen
  Arraylaenge. */
int doubleVec(FILE *fp, double *m) {
  /*Anzahl der enthaltenen Zahlen und Zeilen*/
  int i=0;

  /*Text-Datei wird zunaechst in "Tokens" zerlegt*/
  char ** lp, ** list = split(NULL, "");
  char buf[BIGBUFSIZE];

  setvbuf(fp, NULL/*buf*/, _IOLBF, BIGBUFSIZE);

  if (! list)
    fputs("no room\n", stderr), exit(1);

  while (fgets(buf, sizeof buf, fp))
    {/* 	++ lno; */
      if (buf[strlen(buf)-1] != '\n') {
/* 	fprintf(stderr, "line %u too long\n", lno); */
	fprintf(stderr, "line too long\n");
	exit(1);
      }
      list = split(list, strsave(buf));
    }

  /*Verarbeite Strings zu doubles*/
  for (lp = list, i=0; * lp; ++ lp, i++) {
    m[i] = atof(*lp);
  }

  return i;
}


/* Gibt double-Matrix m in Datei aus*/
void printDoubleMat(FILE *fp, double *A[], int m, int n) {
  int i,j;
  for(i=0; i<m; i++) {
    for(j=0; j<n; j++) {
      fprintf(fp, "%f\t", A[i][j]);
    }
    fprintf(fp, "\n");
  }
}

/* Gibt double-Vektor m in Datei aus*/
void printDoubleVec(FILE *fp, double *v, int n) {
  int i;
  for(i=0; i<n; i++) {
    fprintf(fp, "%f\t", v[i]);
  }
  fprintf(fp, "\n");
}

/*Liest eine Matrix mit int-Werten aus einer Datei ein. Es wird
  zunaechst von einer quadratischen Matrix ausgegangen, Speicher muss
  reserviert sein.*/
void intMat(FILE *fp, int *m[], int M, int N) {
  /*Anzahl der enthaltenen Zahlen und Zeilen*/
  int i,j;
  for (i=0; i< M; i++) {
    for(j=0; j<N; j++) {
      int n;
      if (!fscanf(fp, "%i", &n)) {
	fprintf(stderr, "Error while reading m[%i][%i]\n", i, j);
	exit(1);
      }
      m[i][j]=n;
    }
  }
}

/*Liest eine Matrix mit int-Werten aus einer Datei ein. Es wird
  zunaechst von einer quadratischen Matrix ausgegangen, Speicher muss
  reserviert sein.*/
void doubleMatNew(FILE *fp, double *m[], int M, int N) {
  /*Anzahl der enthaltenen Zahlen und Zeilen*/
  int i,j;
  for (i=0; i< M; i++) {
    for(j=0; j<N; j++) {
      double d;
      if (!fscanf(fp, "%lf", &d)) {
	fprintf(stderr, "Error while reading m[%i][%i]\n", i, j);
	exit(1);
      }
      m[i][j]=d;
    }
  }
}

/* Reads n doubles from string s. Returns number of doubles read. */
int sgetDoubles(double *d, int n, const char *s) {
  int i=0;
  char *cp;
  for (cp = strtok(s, DELIM); cp && (i<n); cp = strtok(NULL, DELIM)) {
    double x;
    /* fprintf(OUT, "sgetDouble(): Token %i: %s\n", i, cp); */
    if(!(sscanf(cp,"%lf", &x))) { 
      fprintf(ERR, "Error reading number %i\n", i);
      return i;
    }
    d[i++]=x;
  }
  return n;
}

/*Liest eine Matrix mit double-Werten aus einer Datei ein. Es wird
  zunaechst von einer MxN-Matrix ausgegangen, Speicher muss reserviert
  sein. Es wird abgebrochen, wenn die Matrix voll ist oder die Datei
  endet. */
int doubleMatVar(FILE *fp, double *m[], int * M, int N) {
  /*Anzahl der enthaltenen Zahlen und Zeilen*/
  int i,j;
  int lc=0;
  printf("%i lines\n", *M);
  for (i=0; i< (*M); i++) {
    /* for(j=0; j<N; j++) { */
    /*   /\* Should read line-wise!!! *\/ */
    /*   double d; */
    /*   int r=fscanf(fp, "%lf", &d); */
    /*   if(r == EOF) {*M=i; return j;} */
    /*   if (!r) { */
    /* 	fprintf(stderr, "Error while reading m[%i][%i]\n", i, j); */
    /* 	exit(1); */
    /*   } */
    /*   m[i][j]=d; */
    /* } */

    char line[LINE_MAX];
    double d;
    if(fgets(line, sizeof line, fp) == NULL) {
      /* fprintf(ERR, "doubleMatVar(): Reached end of file\n"); */
      *M=i;
      /* This means that last line was read successfully until the end */
      return N;
    }
    /* fprintf(OUT, "doubleMatVar(): line %i:\n%s", i, line); */
    if( (j=sgetDoubles(m[i], N, line)) < N) {
      *M=i;
      return j;
    }
    /* for(j=0; j<N; j++) { */
    /* /\* for (cp = strtok(s, DELIM); cp; cp = strtok(NULL, DELIM)) { *\/ */

    /*   /\* } *\/ */
    /*   int e=sscanf(line, "%lf", &d); */
    /*   fprintf(OUT, "doubleMatVar(): sscanf=%i, d=%f\n", e,d); */
    /*   if(!e) { */
    /* 	fprintf(ERR, "doubleMatVar(): Too few values in line %i\n", i); */
    /* 	*M=i; */
    /* 	return j; */
    /*   } */
    /*   m[i][j]=d; */
    /* } */
  }
  return N;
}



/*Liest einen Vektor von N int-Werten aus einer Datei ein. Speicher
  muss reserviert sein, gibt Anzahl tatsaechlich gelesener Werte
  zurueck.*/
int intVec(FILE *fp, int *m, int N) {
  /*Anzahl der enthaltenen Zahlen und Zeilen*/
  int i;
  for (i=0; i< N; i++) {
    int n;
    if (fscanf(fp, "%i", &n)<1) break;/*  { */
/*       fprintf(stderr, "Error while reading m[%i]\n", i); */
/*       exit(1); */
/*     } */
    else  {fprintf(stdout, "%i\n", n); m[i]=n;}
  }
  return i;
}



#define MAXCOL 50
#define MAXROW 6000

/* Wandelt String schrittweise in doubles um, reserviert Speicher fuer
   v, wenn noetig, speichert Laenge von v in nV. */
double* str2doubles(char *s, int *nV) {
  char * cp, *endptr;
  double * lp;
  unsigned lim;
  double *v;

  /* if (v) { */
  /*   /\* lim bis zum Ende von list bewegen, um zu verlaengeren, lim ist */
  /*      die Anzahl der Eintraege *\/ */
  /*   for (lp = v; * lp; ++ lp); */
  /*   lim = lp - v + 1; */
  /* } */
  /* else
 */ 
  if ((v = malloc(MAXCOL*sizeof(double))))
    lp = v, lim = MAXCOL;
  /* kein Speicher vorhanden */
  else return NULL;

  /* cp enthaelt in jedem Schritt einen Zeiger auf das naechste
     Token... und lp enthaelt jeweils das nach double umgewandelte
     cp */
  for (cp = strtok(s, DELIM) ; cp; cp = strtok(NULL, DELIM)) {
    if (lp - v == (lim - 1)) {
      if (! (lp = realloc(v, (lim + MAXCOL)*sizeof(double) ))) {
	/* v[lim-1] = NULL; */
	*nV=lim-1;
	return v;
      }
      v = lp, lp = v + lim - 1, lim += MAXCOL;
    }
    /* Umwandlung */
    *lp=strtod(cp, &endptr);
    lp++;
  }
  /* How many numbers were converted?*/
  *nV=lp-v;
  return realloc(v, (lp - v)*sizeof(double));
}

/* Liest zeilenweise Textdatei unter der Annahme, dass die Anzahl der
   Spalten von der ersten Zeile an gleichbleibt. */
double ** file2mat(FILE *fp, int *nR, int *nC) {
  double **out = malloc(MAXROW*sizeof(double *));
  char buf[BIGBUFSIZE];
  int lim=MAXROW;
  int nCcurr;

  setvbuf(fp, NULL, _IOLBF, BIGBUFSIZE);
  fgets(buf, sizeof buf, fp);
  out[0]=/* dblsave( */
    str2doubles(buf, nC);
		 /* ); */
  *nR=1;
  while (fgets(buf, sizeof buf, fp)) {
    if (buf[strlen(buf)-1] != '\n') {
	fprintf(stderr, "line %u too long\n", *nR++);
	exit(1);
    }
    if(*nR == lim-1) {
      if (! (out = realloc(out, (lim + MAXROW)*sizeof(double *) ))) {
	/* v[lim-1] = NULL; */
	*nR=lim-1;
	return out;
      }
      lim+=MAXROW;
    }
    /* fprintf(stdout, "Line %i\n", *nR); */
    out[*nR] = /* dblsave( */str2doubles(buf, &nCcurr)/* ) */;
    *nR+=1;
  }
  return realloc(out, (*nR)*sizeof(double *));
}

/* Gibt Spalte j einer dbl-Matrix zurueck */
double *dblMatCol(double **m, int j, int nR) {
  int i;
  double *m_J=malloc(nR*sizeof(double));
  for(i=0; i<nR; i++) {
    m_J[i]=m[i][j];
  }
  return m_J;
}

/* Gibt Zeile i einer dbl-Matrix zurueck */
double *dblMatRow(double **m, int i, int nC) {
  int j;
  double *m_I=malloc(nC*sizeof(double));
  /* for(j=0; j<nC; j++) { */
  /*   m_I[i]=m[i][j]; */
  /* } */
  m_I=memcpy(m_I, m[i], nC*sizeof(double));
  return m_I;
}

/* selects nCol columns starting from j0 from matrix A and writes them
   to long vector v */
void selectColumns(const double *A[], 
		   int j0, int nCol,
		   double *v, int M) {
  int i, j,k=0;

  for(j=j0; j<j0+nCol; j++) {
    for(i=0; i<M; i++) {
      v[k++]=A[i][j];
    }
  }
}
