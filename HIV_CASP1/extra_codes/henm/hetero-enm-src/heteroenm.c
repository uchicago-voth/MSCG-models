/*
 * This program is used to do normal mode analysis (NMA) by an elastic
 * network model (ENM). By Zhiyong "John" Zhang, 02/28/2008.
 * This version is cleaned up by AVS in Aug 2011.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "pca.h"

extern void read_x(FILE *file, int natm, rvec *x);
extern int dsyevr_(char *, char *, char *, long int *, double *,
                   long int *, double *, double *, long int *,
                   long int *, double *, long int *, double *,
                   double *, long int *, long int *, double *,
                   long int *, long int *, long int *, long int *);

int rdpara(int argc, char *argv[], FILE **fp)
{
 int i;

 for (i=1; i<argc; i++)
 {
    if (strcmp(argv[i], "-icf") == 0)
    {
      fp[0] = fopen(argv[++i], "r");
    }
    else if (strcmp(argv[i], "-ibn") == 0)
    {
      fp[1] = fopen(argv[++i], "r");
    }
    else if (strcmp(argv[i], "-ipr") == 0)
    {
      fp[2] = fopen(argv[++i], "r");
    }
    else if (strcmp(argv[i], "-orf") == 0)
    {
      fp[3] = fopen(argv[++i], "w");
    }
 }
 return(1);
}

int main(int argc, char *argv[])
{
 FILE *ircf, *idis, *ibon, *ikfc, *ipar, *orsf, *fp[20];
 int i, idum, j, jdum, k, n, m, i2, j2, k2, i3, j3, l2, natm, nmodes, **bon;
 long int ndim, il, iu, lm, lwork, liwork, info, *isuppz, *iwork;
 double fc, vl, vu, abstol, rmsf, sprj, prj, T;
 double *hm, *work, *val, *vec, **dis, **kr, d,ki;
 rvec *x, *xvec;
 char str[STRLEN];

 if (argc == 1 || strcmp(argv[1], "-h") == 0)
 {
   fprintf(stderr, "usage: %s -icf coordinates.xyz -ibn cgk.dat -ipr parameters.in -orf enm-bond-flucs.dat\n", argv[0]);
   exit(1);
 }

 if (rdpara(argc, argv, fp) == 0) exit(0);
 ircf = fp[0];
 ibon = fp[1];
 ipar = fp[2];
 orsf = fp[3];

 /* control parameters */
 fgets2(str, STRLEN, ipar); // one-line statement of the file

 /* number of atoms */
 fgets2(str, STRLEN, ipar); // one-line comment
 fgets2(str, STRLEN, ipar);
 sscanf(str, "%d", &natm);

 /* the last number of normal mode we want to use */
 fgets2(str, STRLEN, ipar); // one-line comment
 fgets2(str, STRLEN, ipar);
 sscanf(str, "%d", &nmodes);

 /* Temperature for kT */
 fgets2(str, STRLEN, ipar); // one-line comment
 fgets2(str, STRLEN, ipar);
 sscanf(str, "%lf", &T);

 ndim = natm*DIM;

 /* read the atomistic coordinates */
 snew(x, natm);
 read_x(ircf, natm, x);

 /* distance matrix */
 snew(dis, natm);
 for (i=0; i<natm; i++)
 {
    snew(dis[i], natm);
 }

 /* read the bond and distance information */
 // initialize arrays
 snew(bon, natm);
 for (i=0; i<natm; i++)
 {
    snew(bon[i], natm);
 }

 for (i=0; i<natm-1; i++)
 {
    for (j=i+1; j<natm; j++)
    {
		bon[i][j]=0;
		bon[j][i]=0;
    }
 }

 /* spring constants */
 snew(kr, natm);
 for (i=0; i<natm; i++)
 {
    snew(kr[i], natm);
 }

 for (i=0; i<natm-1; i++)
 {
    for (j=i+1; j<natm; j++)
    {
		kr[i][j]=0;
		kr[j][i]=0;
    }
 }
 
 /* read the bonds file */
 while(fgets2(str, STRLEN, ibon)!=NULL){
    sscanf(str, "%d%d%lf%lf", &i,&j, &d, &ki);
	i-=1;
	j-=1;
	kr[i][j] = 2*4184/(T*1.38*6.022)*ki;
              // original ki are measured in kcal/mol/Angstrom^2
              // kr[][] are dimensionalless parameters, equal to the physical k_ij divided by (kT/2)
	kr[j][i] = kr[i][j];
	dis[i][j]=d;
	dis[j][i] = dis[i][j]; 
	bon[i][j] = 1;
	bon[j][i] = 1;

 }

 /* Hessian matrix */
 // "hm" here is the Hessian matrix (the matrix of the second derivatives of the energy w.r.t. the coordinates) divided
 // by kT. The elements of hm are measured in Angstrom^(-2).
 snew(hm, ndim*ndim);
 for (i=0; i<ndim; i++)
 {
    for (j=i; j<ndim; j++)
    {
       i2 = i/3;
       j2 = j/3;
       k2 = i%3;
       l2 = j%3;
       if (j2 != i2)
       {
         if( bon[i2][j2]==1)  hm[ndim*j+i] = -kr[i2][j2]*(x[i2][k2]-x[j2][k2])*(x[i2][l2]-x[j2][l2])/(sqr(dis[i2][j2]));
    	 //fprintf(stderr,"%i %i %f\n",i2, j2, hm[ndim*j+i]);
       }
       else
       {
         for (j3=0; j3<natm; j3++)
         {
            if (j3 != i2 && bon[i2][j3]==1)
            {
              hm[ndim*j+i] += kr[i2][j3]*(x[i2][k2]-x[j3][k2])*(x[i2][l2]-x[j3][l2])/(sqr(dis[i2][j3]));
              //fprintf(stderr,"%i %i %f\n",i2,j2, hm[ndim*j+i]);
	    }
         }
       }
       hm[ndim*i+j] = hm[ndim*j+i]; // symmetric Hessian matrix
    }
 }

 /* free some space */
 for (i=0; i<natm; i++)
 {
    sfree(kr[i]);
 }

 /* diagonalization */
 /* allocate space */
 vl = 0.0;
 vu = 1000.0;
 il = 7; // ignore the first six zero-eigenvalue modes
 iu = nmodes;
 abstol = 0.0;
 lwork = 30*ndim;
 liwork = 15*ndim;
 snew(val, ndim);
 snew(vec, ndim*nmodes);
 snew(isuppz, 2*nmodes);
 snew(work, lwork);
 snew(iwork, liwork);

 dsyevr_("V", "I", "U", &ndim, hm, &ndim, &vl, &vu, &il, &iu,
         &abstol, &lm, val, vec, &ndim, isuppz, work, &lwork, iwork,
         &liwork, &info);

 snew(xvec, natm);

 /* calcualate ENM fluctuations onto the bonds */
 // the formula used here come from the lowest-order expansion of |r_ij| around its equilibrium value
 for (i=0; i<natm-1; i++)
 {
    for (j=i+1; j<natm; j++)
    {
       prj = 0.0;
       if (bon[i][j] == 1)
       {
         for (n=0; n<nmodes-6; n++)
         {
            sprj = 0.0;
            for (m=0; m<DIM; m++)
            {
               xvec[i][m] = vec[ndim*n+i*DIM+m];
               xvec[j][m] = vec[ndim*n+j*DIM+m];
               //sprj += ((x[i][m]-x[j][m])/dis[i][j])*(xvec[i][m]/sqrt(val[n])-xvec[j][m]/sqrt(val[n]));
               /* In theory, val[n] must be positive. But sometimes it turns out to be -1.e-8 or like that due
                  to machine error. This must be corrected */
               if (val[n]<-1.e-6)
               {
                 fprintf(stderr, "Warning from hetero-enm: eigenvalue for n=%d was %f, now replaced by +1.e-8\n", n, val[n]);
               }
               if (val[n]<=0) val[n]=1.e-8;
               sprj += ((x[i][m]-x[j][m])/dis[i][j])*(xvec[i][m]-xvec[j][m])/sqrt(val[n]);
            }
            prj += sqr(sprj);
         }
         fprintf(orsf, "%4d %4d %15.9lf\n", i+1, j+1, sqrt(prj));
       }      
    }
 }

 /* close file */
 fclose(ircf);
 fclose(ibon);
 fclose(ipar);
 fclose(orsf);

 /* the end of main */
 return 0;
}
