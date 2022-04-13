#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "rk4.h"
#include <errno.h>

#define N 300 //r grid
#define NZ 200 //z grid
#define DL 0.631 //d-m interaction coef
#define DP 0.005//dumping coef
#define AN 1  //anisotropy coef
#define V -0.245 // velocity
#define PHI M_PI // stereographic projection angle
#define T_STOP 35 // runge kutta final time
#define A 4 // radius
#define INIT 0 //

int magnIndex(int i, int j);
double deriv2r(double *m,int i,double dr);
double deriv1r(double *m,int i, double dr);
double  deriv2z(double *m,int i,double dz);
double deriv1r(double *m,int i,double dz);
void laplace(double *m,double *r,double *z,double *laplacian,int ri,int zi,double dr,double dz);
void vortexRing_init_values(double *m0,double *r, double *z,double a,double phi);
void LLequation(double t, double *m,double *dm,double *r,double *z,double dr,double dz);
void boundary_conditions(double *m);
double EXenergy(double *m, double *r,double *z, double dr,double dz);
double ANenergy(double *m, double *r,double *z, double dr,double dz);
double DMenergy(double *m, double *r,double *z, double dr,double dz);
void DMenergy_Density(double *m, double *r,double *z, double dr,double dz,double a,double phi);
void center_angle(double *y,double *r,double *z, double dr, double dz);
void DMenergy_zDensity(double *m, double *r,double *z, double dr,double dz,double a,double phi);
void skyrmion_init_values(double *m0,double *r,double *z);
void boundary_conditions(double *m);
void moving_LLequation(double t, double *m,double *dm,double *r,double *z,double dr, double dz);
double mz_integral(double *m,double *r,double *z, double dr, double dz);
