#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "pca.h"

void *save_calloc(char *name,unsigned nelem,unsigned elsize)
{
  void *p;
                                                                                
  p = NULL;
  if ((nelem == 0) || (elsize == 0))
    p = NULL;
  else
    p = calloc((size_t)nelem,(size_t)elsize);
  return p;
}
                                                                                
void save_free(char *name,void *ptr)
{
  if (ptr != NULL)
    free(ptr);
}

char *fgets2(char *line, int n, FILE *stream)
/* This routine reads a string from stream of max length n
 * and zero terminated, without newlines
 * line should be long enough (>= n)
 */
{
  char *c;
  if (fgets(line,n,stream)==NULL) return NULL;
  if ((c=strchr(line,'\n'))!=NULL) *c=0;
  return line;
}

void read_x(FILE *file, int natm, rvec *x)
{
 int i, j;
 char str[STRLEN];
 double d[DIM];
 char *b;
 snew(b,10);
 // read two header lines
 fgets2(str, STRLEN, file);
 fgets2(str, STRLEN, file);                                                                               
 for (i=0; i<natm; i++) {
    fgets2(str, STRLEN, file);
    sscanf(str, "%s%lf%lf%lf", b, &d[0], &d[1], &d[2]);
    
    for (j=0; j<DIM; j++)
       x[i][j] = d[j];
 }
}
