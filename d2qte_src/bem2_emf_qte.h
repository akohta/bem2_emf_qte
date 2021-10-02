#ifndef BEM2_EMF_QTE_H_
#define BEM2_EMF_QTE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <mkl.h>

#include "physical_constant.h"
#include "osu_mksa.h"
#include "my_utils.h"
#include "gauleg.h"
#include "incident_field.h"
#include "pb_elem.h"


// -- d2qte_setup.c --
void read_data(int argc,char **argv,DOMD *md);    // read datafile
void print_data(DOMD *md);                        // print loaded data
void print_data_MKSA(DOMD *md);                   // print loaded data in MKSA system of units
void initialize_domd(DOMD *md);                   // memory allocation and initialize the data
void mfree_domd(DOMD *md);                        // memory free
int domain_id(double *rt,DOMD *md);               // return domain id at point rt, rt=(x,y)
void dat_read (char *dname,DOMD *md);             // read datafile output by dat_write()
void dat_write(char *dname,DOMD *md);             // write datafile
void output_node_particles(char *fname,DOMD *md); // outputs the nodes as point cloud data ( .particles file ) 


// -- d2qte_solve_bieq.c --
void solve_bieq(DOMD *md);                        // solve boundary integral equations, for bv_solver


// -- d2qte_field.c --
int Ez_s(double complex *Ez,double *rt,int type,DOMD *md); // scattered or internal field
int Ez_t(double complex *Ez,double *rt,int type,DOMD *md); // total field 
int Ez_i(double complex *Ez,double *rt,int type,DOMD *md); // incident field
int Hr_s(double complex *Hr,double *rt,int type,DOMD *md); // scattered or internal field
int Hr_t(double complex *Hr,double *rt,int type,DOMD *md); // total field
int Hr_i(double complex *Hr,double *rt,int type,DOMD *md); // incident field
int EH_s(double complex *Ez,double complex *Hr,double *rt,int type,DOMD *md); // scattered or internal field
int EH_t(double complex *Ez,double complex *Hr,double *rt,int type,DOMD *md); // total field
int EH_i(double complex *Ez,double complex *Hr,double *rt,int type,DOMD *md); // incident field
// optputs
// Ez=E_z, Hr[0]=H_x,Hr[1]=H_y, return domain id. 
// inputs 
// rt[0]=x, rt[1]=y, type : select integration method, md : pointer of DOMD object.
// type=0 : 3-point Gauss-Legendre, type=1 : 7 point Guass-Kronrod (Gauss-Kronrod extension of 3-point GL)
// type=2 : GLH (defined in d2te_const.h) point Gauss-Legendre, type=3 : DE method (for test)
// type>3 : integration with numerical validation 

// boundary value 
double complex Ez_bv(int did,double eta_t,int t,DOMD *md);
void EH_bv(double complex *Ez,double complex *Hr,int did,double eta_t,int t,DOMD *md);
// inputs
// did : domain id, eta_t : parameter of point on each element ( -1 ~ 1 ), t : element id, pointer of DOMD object.

// tangential directinal derivative of boundary vale at the node
double complex dEzdt_bv_node_ndmtd(int did,int t,int tn,DOMD *md); // using numerical differentiation
double complex dEzdt_bv_node_dbieq(int did,int t,int tn,DOMD *md); // using derivative boundary integral equation
// return tangential directinal derivative of Ez at the node.
// inputs
// did : domain id, t : element id, tn : node id (0~2), md : pointer of DOMD object.


// -- d2qte_force.c --
void force_FN(double *Fr,double *Nz,double *rc,int type,DOMD *md);
// outputs
// Fr[0]=F_x, Fr[1]=F_y (radiation force), Nx=N_z (radiation torque)
// inputs
// rc[0]=x,rc[1]=y : coordinate of rotation center,
// type : select integration medthod, md : pointer of DOMD object.
// type=0 : 3-point Gauss-Legendre, type>0 : 7-point Gauss-Kronrod.

#endif
