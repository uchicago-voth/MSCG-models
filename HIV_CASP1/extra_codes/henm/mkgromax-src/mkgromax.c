/*
 * This program prepares input files for Gromacs.
 * Written by Andrea Grafmueller.
 * This version is cleaned up in Aug 2011, AVS
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "pca.h"


typedef int bool; 
#define true (1)
#define false (0)


extern void read_x(FILE *file, int natm, rvec *x);

void pseudoinverse(double *l, double  **v, double **crc, int n, int nm) ;

int rdpara(int argc, char *argv[], FILE **fp)
{
 int i;

 for (i=1; i<argc; i++)
 {
    if (strcmp(argv[i], "-ipr") == 0)
    {
      fp[0] = fopen(argv[++i], "r");
    }
    else if (strcmp(argv[i], "-imw") == 0)
    {
      fp[1] = fopen(argv[++i], "r");
    }
    else if (strcmp(argv[i], "-iat") == 0)
    {
      fp[2] = fopen(argv[++i], "r");
    }
    else if (strcmp(argv[i], "-ibo") == 0)
    {
      fp[3] = fopen(argv[++i], "r");
    }
    else if (strcmp(argv[i], "-opr") == 0)
    {
      fp[4] = fopen(argv[++i], "w");
    }
    else if (strcmp(argv[i], "-otp") == 0)
    {
      fp[5] = fopen(argv[++i], "w");
    }
    else if (strcmp(argv[i], "-ocr") == 0)
    {
      fp[6] = fopen(argv[++i], "w");
    }
    else if (strcmp(argv[i], "-oat") == 0)
    {
      fp[7] = fopen(argv[++i], "w");
    }
 }
 fp[8] = fopen("ffcharmm27nb.itp","w");
 return(1);
}


int main(int argc, char *argv[])
{
 FILE *icor, *imas, *ipar, *ibonds, *oitp, *otop, *ocrd, *oatp,*offnb, *fp[20];
 int i, j, k, m, nt, i2, j2, k2, i3, j3, l2, natm, nbonds, nAtomTypes;
 long int ndim, *vBondEnds;
 double *mass, *mtypes;
 double *vBondLK;
 bool **bondij, add=false;
 rvec *x;
 char str[STRLEN], **name, atname[5], **types, *itpname, mypath[STRLEN];


 // for now set itpname here
 itpname = "enm.itp";

 if (argc == 1 || strcmp(argv[1], "-h") == 0)
 {
   fprintf(stderr, "usage: %s  -imw mass.dat -ipr gro.in -iat cg.xyz -ibo bond.dat -opr enm.itp -otp enm.top -ocr enm.gro -oat enm.atp \n", argv[0]);
   exit(1);
 }

 if (rdpara(argc, argv, fp) == 0) exit(0);
 
 ipar = fp[0];
 imas = fp[1];
 icor = fp[2];
 ibonds = fp[3];
 oitp = fp[4];
 otop = fp[5];
 ocrd = fp[6];
 oatp = fp[7]; 
 offnb=fp[8];

 /* control parameters */
 fgets2(str, STRLEN, ipar); // one-line statement of the file

 /* number of sites */
 fgets2(str, STRLEN, ipar); // one-line comment
 fgets2(str, STRLEN, ipar);
 sscanf(str, "%d", &natm);

 /* number of sites */
 fgets2(str, STRLEN, ipar); // one-line comment
 fgets2(str, STRLEN, ipar);
 sscanf(str, "%d", &nbonds);

 // number of different site types (get fromm script)
 fgets2(str, STRLEN, ipar); // one-line comment
 fgets2(str, STRLEN, ipar);
 sscanf(str, "%d", &nAtomTypes);

 //the include path
 fgets2(str, STRLEN, ipar); // one-line comment
 fgets2(str, STRLEN, ipar);
 sscanf(str, "%s", mypath);

 ndim = natm*DIM;
 
 /* read the coordinates */
 snew(x, natm);
 read_x(icor, natm, x);
 //    fprintf(stderr, "%lf  %lf  %lf" , x[1][0], x[1][1], x[1][0]);


 /* read the residue masses */
 snew(mass, natm);
 snew(name, natm);
 snew(types, nAtomTypes);
 snew(mtypes, nAtomTypes);
 for(i=0; i<natm; i++){
	snew(name[i],10);
 }
 for(i=0; i<nAtomTypes; i++){
 	snew(types[i],10);
 }

 nt=0;

 for (i=0; i<natm; i++)
 {
    fgets2(str, STRLEN, imas);
    sscanf(str, "%s%lf", name[i], &mass[i]);
    add = true;
    for(j=0; j<nt; j++){
	if (strcmp(name[i],types[j]) == 0) add=false;
    }    
    if(add){
	 strcpy(types[nt],name[i]);
	 mtypes[nt] = mass[i];
    }
    nt ++;
    if (nt > nAtomTypes) fprintf(stderr, "error: more atom types than expected\n");
 }

 /* read the bondlengths and force constants */
 snew(vBondEnds, 2*nbonds);
 snew(vBondLK, 2*nbonds);
 snew(bondij, natm);
 for(i=0; i< natm; i++){
   snew(bondij[i], natm);
 }
 for(i=0; i<natm; i++) for(j=i; j<natm; j++)
 {
	bondij[i][j] = false;
	bondij[j][i] = false;
 }

 for (i=0; i<nbonds; i++)
 {
    fgets2(str, STRLEN, ibonds);
    sscanf(str, "%li%li%lf%lf", &vBondEnds[2*i], &vBondEnds[2*i+1], &vBondLK[2*i], &vBondLK[2*i+1]);
    bondij[vBondEnds[2*i]-1][vBondEnds[2*i+1]-1]=true;
    bondij[vBondEnds[2*i+1]-1][vBondEnds[2*i]-1]=true;
 }
 

 //**** done reading input ************************************************************* 
 
 //print the gromacs input

 //atp file 
 for(i=0; i< nAtomTypes; i++)
 {
    fprintf(oatp,"%s  %5.4f\n",types[i], mtypes[i] );
 }

 // itp files
 fprintf(offnb,"[ atomtypes ]\n;  name  at. num  mass charge  ptype  sig  eps\n");
 for (i=0; i<nAtomTypes; i++)
 {
    fprintf(offnb,"%4s%6i%10.4f%8.2f%4s%8.2f%8.2f\n",name[i],i+1,mass[i],0.0,"A",0.0,0.0);
 } 
 fprintf(oitp, ";\t itp file for enm calculation \n\n");
 fprintf(oitp,"[ moleculetype ]\n");
 fprintf(oitp,";  Name    nreexcl\n");
 fprintf(oitp,"ENM_mol       0\n\n");
 fprintf(oitp,"[ atoms ]\n;  nr  type  resnr  residue atom  cgnr  charge  mass\n");
 for (i=0; i< natm; i++){
    fprintf(oitp,"%4i%6s%6i%6s   C%-5i%5i%8.3f  %8.4f\n",i+1,name[i], 1,"ENM", i+1,i+1,0.0, mass[i]);
 } 
 fprintf(oitp,"\n[ bonds ]\n; ai   aj   func  l0  c1\n");
 for(i=0; i<nbonds; i++){
    i2 = vBondEnds[2*i]-1;     
    j2 = vBondEnds[2*i+1]-1;	
    fprintf(oitp,"%5i%5i%6i%12.6f%12.6f\n",i2+1,j2+1,1, vBondLK[2*i]/10, 836.8*vBondLK[2*i+1]);    
 }
 
 //topology file 
 fprintf(otop, ";\t topology file for enm calculation \n\n;  Include forcefield parameters\n\n");
 fprintf(otop, "#include  \"%s/ffcharmm27.itp\"\n\n",mypath);
 fprintf(otop, ";  Include Molecule topology\n\n");
 fprintf(otop, "#include \"%s\"\n\n\n", itpname);
 fprintf(otop, "[ System ]\n; Name \n ENMsystem\n\n");
 fprintf(otop, "[ molecules ]\n; Compound       #mols\n");
 fprintf(otop, " ENM_mol 	1");
 
 // gro - coordinate file
 
 fprintf(ocrd, "Funny Gromacs Commentary\n");
 fprintf(ocrd,"%10i\n",natm);
 for (i=0; i< natm; i++){
    sprintf(atname,"C%i",i+1);
    fprintf(ocrd,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",1,"ENM",atname,i+1,x[i][0]/10,x[i][1]/10,x[i][2]/10);
 }
 fprintf(ocrd, " 80  80  80");
 

 /* close files */
 fclose(ipar);
 fclose(imas);
 fclose(icor);
 fclose(ibonds);
 fclose(oatp);
 fclose(oitp);
 fclose(otop);
 fclose(ocrd);

 /* the end of main */
 return 0;
}
