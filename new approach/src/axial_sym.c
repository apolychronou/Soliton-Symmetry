#include "axialSym.h"

int main(){
  // double *y;
  // double m0[3*N*NZ]={0};
  // double *m0;
  int i=0;
  double dr=0.1;
  double dz=0.1;
  double r[N];
  double z[NZ];
  int steps;
  double *t;
  double tspan[2];
  double a=A,phi=PHI;
  double *y;
  // double *m,*prevm;

  for (i=0;i<N;i++){
    r[i]=i*dr;
  }
  for (i=0;i<NZ;i++){
    z[i]=-NZ/2*dz+i*dz;
  }

  tspan[0]=0.0;
  tspan[1]=T_STOP;

  steps=(tspan[1]-tspan[0])/(dr*dr*0.1);

  t = ( double * ) malloc ( ( steps + 1 ) * sizeof ( double ) );
  y=(double*)malloc((3*N*NZ*2)*sizeof(double));
  // m=(double*)malloc((3*N*NZ)*sizeof(double));
  // prevm=(double*)malloc((3*N*NZ)*sizeof(double));
  for (i=0;i<3*N*NZ*2;i++){
    y[i]=0;
  }


  if(INIT){
    center_angle(y,r,z,dr,dz);
    DMenergy_Density(y,r,z,dr,dz,a,phi);
    DMenergy_zDensity(y,r,z,dr,dz,a,phi);
    return 0;
  }

  vortexRing_init_values(y,r,z,a,phi);
  // skyrmion_init_values(y,r,z);
  // m0=y;

  clock_t begin = clock();
  // rk4 ( moving_LLequation, tspan, y, steps, 3*N*NZ, t, y,r,z,dr,dz );
  rk4 ( LLequation, tspan, y, steps, 3*N*NZ, t, y,r,z,dr,dz );

  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  printf("execution time: %lf\n",time_spent);

  free(t);
  free(y);
  // free(m);
  // free(prevm);
  return 0;
}
