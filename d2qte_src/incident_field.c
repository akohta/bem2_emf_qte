/*
 * incident_field.c
 *
 *  Created on: Sep 12, 2018
 *      Author: ohta
 */
 
#include "incident_field.h"


void read_infd(char *fname,INFD *wd)
{
  FILE *fp;
  if((fp=fopen(fname,"rt"))==NULL){    printf("Can not open the %s file. Exit...\n",fname);    exit(1);  }
  strcpy(wd->fname,fname);
  char buf[256]="";   double tmpd,tmpd2;
  if(fgets(buf,256,fp)==NULL){
    printf("incident_field.c, read_infd(), failed to read the line. exit...\n");
    exit(1);
  }
  if(fgets(buf,256,fp)==NULL){
    printf("incident_field.c, read_infd(), failed to read the line. exit...\n");
    exit(1);
  }

  if(fscanf(fp,"%lf",&tmpd)!=1){
    printf("incident_field.c, read_infd(), failed to read the real(E0). exit...\n");
    exit(1);
  }
  if(fscanf(fp,"%lf",&tmpd2)!=1){
    printf("incident_field.c, read_infd(), failed to read the imag(E0). exit...\n");
    exit(1);
  }
  wd->E0      =tmpd+I*tmpd2;
  if(fscanf(fp,"%lf",&tmpd)!=1){
    printf("incident_field.c, read_infd(), failed to read the lambda0. exit...\n");
    exit(1);
  }    
  wd->lambda0 =tmpd;
  if(fscanf(fp,"%lf",&tmpd)!=1){
    printf("incident_field.c, read_infd(), failed to read the ne. exit...\n");
    exit(1);
  }
  wd->ne      =tmpd;
  if(fscanf(fp,"%lf",&tmpd)!=1){
    printf("incident_field.c, read_infd(), failed to read the angle. exit...\n");
    exit(1);
  }
  wd->angle   =tmpd;
  fclose(fp);

  // init data
  wd->k0=2.0*M_PI/wd->lambda0;
  wd->kxh=cos(wd->angle);
  wd->kyh=sin(wd->angle);
}

void print_infd(INFD *wd)
{
  printf("-- incidnt field data --\n");
  printf("filename                 : %s\n",wd->fname);
  printf("planewave amplitude E0   :%7.6g+%7.6gI\n",creal(wd->E0),cimag(wd->E0));
  printf("wavelength in vacuum     : %15.14g\n",wd->lambda0);
  printf("refractive index         : %15.14g\n",wd->ne);
  printf("incident angle     [rad] : %15.14g\n",wd->angle);
}

void print_infd_MKSA(INFD *wd)
{
  printf("-- incidnt field data MKSA--\n");
  printf("filename                   : %s\n",wd->fname);
  printf("planewave amplitude E0[V/m]:%7.6g+%7.6gI\n",creal(OSUtoMKSA_ElectricField(wd->E0)),cimag(OSUtoMKSA_ElectricField(wd->E0)));
  printf("wavelength in vacuum    [m]: %15.14g\n",OSUtoMKSA_length(wd->lambda0));
  printf("refractive index           : %15.14g\n",wd->ne);
  printf("incident angle       [rad] : %15.14g\n",wd->angle);
}

void infd_grad(double complex *Ez,double complex *gradEz,double *r,INFD *wd)
{
  double theta=wd->k0*wd->ne*(wd->kxh*r[0]+wd->kyh*r[1]);
  *Ez=wd->E0*(cos(theta)+I*sin(theta));
  gradEz[0]=I*wd->k0*wd->ne*wd->kxh*(*Ez);
  gradEz[1]=I*wd->k0*wd->ne*wd->kyh*(*Ez);
}

double complex infd_Ez(double *r,INFD *wd)
{
  double theta=wd->k0*wd->ne*(wd->kxh*r[0]+wd->kyh*r[1]);
  return wd->E0*(cos(theta)+I*sin(theta));
}

void infd_Hr(double complex *Hr,double *r,INFD *wd)
{
  double theta=wd->k0*wd->ne*(wd->kxh*r[0]+wd->kyh*r[1]);
  double complex ch=wd->E0*wd->ne;
  double complex ex=cos(theta)+I*sin(theta);
  Hr[0]= ch*wd->kyh*ex;
  Hr[1]=-ch*wd->kxh*ex;
}

void infd_EH(double complex *Ez,double complex *Hr,double *r,INFD *wd)
{
  double theta=wd->k0*wd->ne*(wd->kxh*r[0]+wd->kyh*r[1]);
  double complex ch=wd->E0*wd->ne;
  double complex ex=cos(theta)+I*sin(theta);
  *Ez=wd->E0*ex;
  Hr[0]= ch*wd->kyh*ex;
  Hr[1]=-ch*wd->kxh*ex;
}
