/*
 * incident_field.h
 *
 *  Created on: Sep 12, 2018
 *      Author: ohta
 */

#ifndef INCIDENT_FIELD_H_
#define INCIDENT_FIELD_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "osu_mksa.h"


typedef struct incident_field_data{
  char fname[128];   // datafile name
  double complex E0; // amplitude
  double lambda0;    // wavelength in vacuum
  double ne;         // refractive index
  double angle;      // incident angle

  double k0;         // wave number in vacuum (=angular frequency in this system)
  double kxh,kyh;    // wave number unit vector
}INFD;

// -- incident_field.c (TE mode plane wave) --
void read_infd(char *fname,INFD *wd); // read incident field data
void print_infd(INFD *wd);            // print incident field data
void print_infd_MKSA(INFD *wd);       // print incident field data in MKSA system of units

void infd_grad(double complex *Ez,double complex *gradEz,double *r,INFD *wd);
// outputs
// Ez=E_z, gradEz[0]=dE_z/dx, gradEz[1]=dE_z/dy
// inputs 
// r[0]=x, r[1]=y,pointer of INFD 
double complex infd_Ez(double *r,INFD *wd); 
// return E_z
void infd_Hr(double complex *H,double *r,INFD *wd); 
// outputs
// H[0]=H_x, H[1]=Hy
void infd_EH(double complex *Ez,double complex *H,double *r,INFD *wd);
// outputs
// Ez=E_z, H[0]=H_x, H[1]=H_y


#endif /* INCIDENT_FIELD_H_ */
