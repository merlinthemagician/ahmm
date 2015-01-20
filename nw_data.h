#ifndef NW_DATA
#define NW_DATA
/*
 * nw_data.h
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

/* Data type for a tokenized line */
typedef struct {
  /* Number of tokens */
  size_t n;
  /* Tokenized lines */
  char **s;
}nw_tokline;

typedef struct {
  /* Number of tokenized lines */
  size_t m;
  /* Array of tokenized lines */
  nw_tokline **lines;
} nw_tokmatrix;

/* Count the number of entries in tl */
int nw_data_nCol(const nw_tokline *tl);

/* Converts string s to tokenizes line nw_tokline. Return NULL if this failed. */
nw_tokline * nw_data_strToTokline(const char *s);

/* Reads file fp line by line and save tokenized lines in tl. Returns number of lines successfully read. */
int nw_data_fpToTokline(FILE *fp, nw_tokline **tl);

/* Print tokenized line */
void nw_data_printTokline(FILE *fp, const nw_tokline *tl);

/* Print tokenized lines */
void nw_data_printTokmatrix(FILE *fp, const nw_tokline **tm, int n);

/* Return N strings in i-th column, starting from row j0 */
char **nw_data_getColumn(const nw_tokline **tm, int i, int j0, int N);

/* SHOULD return index until which conversion was successful */
/* Convert N strings in i-th column, starting from row j0 to doubles */
int nw_data_columnToDouble(const nw_tokline **tm, int i, int j0,  int N,
			   double *v);

/* Converts columns nCol columns starting from index col0 to
   double. In each column the first entry is j0.
*/
int nw_data_convertColumnsToDouble(const nw_tokline **tm, int col0, int nCol,
				   int j0,  int N, double *v);

double *nw_data_selectDoublesByRowColumns(const nw_tokline **tm, 
					  const char *rh, const char *ch[], 
					  int nCols, int N);
#endif
