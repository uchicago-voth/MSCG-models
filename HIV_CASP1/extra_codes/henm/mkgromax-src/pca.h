/* Some parameters and functions needed by LFA, from GROMACS.
 * By Zhiyong Zhang, 02/02/2004.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "vec.h"

#define XX 0
#define YY 1
#define ZZ 2
#define STRLEN 200
#define snew(ptr,nelem) (ptr)=save_calloc(#ptr,(nelem),sizeof(*(ptr)))
#define sfree(ptr) save_free(#ptr,(ptr))

void *save_calloc(char *name,unsigned nelem,unsigned elsize);
void save_free(char *name,void *ptr);

extern char *fgets2(char *s, int n, FILE *stream);

static inline void oprod(rvec a,rvec b,rvec c)
{
  c[XX]=a[YY]*b[ZZ]-a[ZZ]*b[YY];
  c[YY]=a[ZZ]*b[XX]-a[XX]*b[ZZ];
  c[ZZ]=a[XX]*b[YY]-a[YY]*b[XX];
}

static inline void clear_mat(matrix a)
{
  a[XX][XX]=a[XX][YY]=a[XX][ZZ]=0.0;
  a[YY][XX]=a[YY][YY]=a[YY][ZZ]=0.0;
  a[ZZ][XX]=a[ZZ][YY]=a[ZZ][ZZ]=0.0;
}
