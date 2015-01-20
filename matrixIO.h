#ifndef MATRIXIO
#define MATRIXIO
/*
 *  matrixIO.h
 *
 * Ivo Siekmann, 30.10.2006
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


/*
 * Ausgabe eines Zeitschritts nach f
 */
void
output(FILE* f, double t, gsl_vector* y);

/*Ausgabe eines gsl_vector*/
void vprint(FILE * stream, gsl_vector * v);

/*Ausgabe eines gsl_vector*/
void vprintfmt(FILE * stream, gsl_vector * v, const char* fmt);

/*Ausgabe eines Zeitschritts*/
void printStep(FILE * stream, double t, double hx, gsl_vector * v);

/* Ausgabe des Zeitschritts i
 * in Datei mit Namen "pre"+Nummerierung
 * digits darf nicht groesser als 9 sein!
 */ 
void matToFile(const gsl_matrix *mat, const char* pre, int i, int digits);

/*Ausgabe einer gsl_matrix in das File f*/
void printMat(FILE *f, const gsl_matrix *m);

/* Ausgabe einer gsl_matrix in das File f. Zeilentrenner rd und
   Spaltentrenner cd sowie Formatstring fmt */
void printMatf(FILE *f, const gsl_matrix *m,
	       char *rd,char *cd, char *fmt);

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
	      int dimx, int dimy);

/*Gibt Loesung als Grauwerte in PGM-Format aus*/
void matToPGM(char *fn, gsl_vector_view *uView,
	      int ind, double s, int dimx, int dimy);

/*Gibt Loesung als double-Werte aus*/
void matToDat(char *fn, gsl_vector_view *uView,
	      int ind, int dimx, int dimy);

/*Ausgabe einer gsl_matrix in das File f*/
void printMatPGM(FILE *f, const gsl_matrix *m, double s);

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
		 double sRot, double sGruen, double sBlau);


/*Ausgabe einer gsl_matrix in das File f*/
void printMatPGM(FILE *f, const gsl_matrix *m, double s);

/*Liest eine Matrix mit double-Werten aus einer Datei ein*/
gsl_matrix *readMat(FILE *fp);

/*Liest einen Vektor mit double-Werten aus einer beliebig
  strukturierten Datei ein. Speicher muss reserviert sein. KEINE Warnung bei
  Ueberschreiten der maximalen Arraylaenge. */
int doubleVec(FILE *fp, double *m);

/*Liest eine Matrix mit double-Werten aus einer Datei ein. Es wird
 *zunaechst von einer quadratischen Matrix ausgegangen.
 */
void doubleMat(FILE *fp, double *m[], int lno);

/* Gibt double-Matrix m in Datei aus*/
void printDoubleMat(FILE *fp, double *A[], int m, int n);

/*Liest eine Matrix mit int-Werten aus einer Datei ein. Es wird von
  einer Matrix mit M Zeilen und N Spalten ausgegangen, Speicher muss
  reserviert sein.*/
void intMat(FILE *fp, int *m[], int M, int n);

/*Liest eine Matrix mit double-Werten aus einer Datei ein. Es wird von
  einer Matrix mit M Zeilen und N Spalten ausgegangen, Speicher muss
  reserviert sein.*/
void doubleMatNew(FILE *fp, double *m[], int M, int n);

/*Liest eine Matrix mit double-Werten aus einer Datei ein. Es wird
  zunaechst von einer MxN-Matrix ausgegangen, Speicher muss reserviert
  sein. Es wird abgebrochen, wenn die Matrix voll ist oder die Datei
  endet. */
int doubleMatVar(FILE *fp, double *m[], int * M, int n);

/* Gibt double-Vektor m in Datei aus*/
void printDoubleVec(FILE *fp, double *v, int n);

/* Counts fields in fp delimited by DELIM */
int countFields(FILE *fp);

/*Liest einen Vektor von N int-Werten aus einer Datei ein. Speicher
  muss reserviert sein.*/
int intVec(FILE *fp, int *m, int n);

/* Wandelt String schrittweise in doubles um, reserviert Speicher fuer
   v, wenn noetig, speichert Laenge von v in nV. */
double* str2doubles(char *s, int *nV);

/* Liest zeilenweise Textdatei unter der Annahme, dass die Anzahl der
   Spalten von der ersten Zeile an gleichbleibt. */
double ** file2mat(FILE *fp, int *nR, int *nC);

/* Liest zeilenweise Textdatei unter der Annahme, dass die Anzahl der
   Spalten von der ersten Zeile an gleichbleibt, speichert den
   gesamten Inhalt in einem Vector.*/
double * file2vec(FILE *fp, int *nR, int *nC);

/* Gibt Zeile i einer dbl-Matrix zurueck */
double *dblMatRow(double **m, int i, int nC);

/* Gibt Spalte j einer dbl-Matrix zurueck */
double *dblMatCol(double **m, int j, int nR);

/* selects nCol columns starting from j0 from matrix A and writes them
   to long vector v */
void selectColumns(const double *A[], 
		   int j0, int nCol,
		   double *v, int M);
#endif
