#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "rk4.h"
#define N 200
#define DL 0.6 //d-m interaction coef
#define DP 0.8 //dumping coef
#define INITIAL 1


float second_derivative(double *m,int i,double dr);
float first_derivative(double *m,int i, double dr);
void laplace(double *m,double *r,double *laplacian,int i,double dr);
void skyrmion_init_values(double *m0,double *r);
void LLequation(double t, double *m,double *dm,double *r,double dr);
double EXenergy(double *m, double *r, double dr);
double ANenergy(double *m, double *r, double dr);
double DMenergy(double *m, double *r, double dr);



int main(int argc, char *argv[]){
  double *y;
  double m0[3*N]={0};
  int i=0;
  double dr=0.1;
  double r[N];
  int steps;
  double *t;
  double tspan[2];
  double m[3*N]={0};
  // double norm;
  double prevm[3*N];

  for (i=0;i<N;i++){
    r[i]=i*dr;
  }
  // for (i=0;i<N;i++){
  //   printf("r:%f",r[i]);
  // }


  tspan[0]=0.0;
  tspan[1]=100.0;
  steps=(tspan[1]-tspan[0])/(dr*dr*0.1);

  t = ( double * ) malloc ( ( steps + 1 ) * sizeof ( double ) );
  y = ( double * ) malloc ( ( steps + 1 ) * 3*N * sizeof ( double ) );

  for(i=0;i<(steps+1)*3*N;i++){
    y[i]=0;
  }
  // laplace(m,r,laplacian,i,dr);

  skyrmion_init_values(m0,r);
  if(INITIAL){
    printf("DM Energy: %lf \n",DMenergy(m0,r,dr));
    for (i=0;i<3*N;i++){
      printf("%lf ",m0[i]);
    }
    return 1;
  }
  clock_t begin = clock();
  rk4 ( LLequation, tspan, m0, steps, 3*N, t, y,r,dr );
  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  for (i=0;i<3*N;i++){
    m[i]=y[steps*3*N+i];
    // printf("m=%f\n",m[i]);
  }

for (i=0;i<3*N;i++){
    prevm[i]=y[(steps-1)*3*N+i];
}

  // for (i=0;i<3*N;i=i+3){
  //   norm=pow(m[i],2)+pow(m[i+1],2)+pow(m[i+2],2);
  //   norm=sqrt(norm);
  //   printf("m:%f\n",norm);
  // }
  double ANener = ANenergy(m,r,dr);
  double DMener = DMenergy(m,r,dr);
  double EXener = EXenergy(m,r,dr);
  double sum= (DMenergy(m,r,dr)+2*ANenergy(m,r,dr));
  printf("Energy: %.10lf ", ANener+DMener+EXener);
  // printf("Prev energy: %.10lf\t",ANenergy(prevm,r,dr)+ DMenergy(prevm,r,dr) + EXenergy(prevm,r,dr));
  printf("Initial energy: %.10lf ",ANenergy(m0,r,dr)+ DMenergy(m0,r,dr) + EXenergy(m0,r,dr));
  printf("sum1= %.10lf ",sum/ANenergy(m,r,dr));
  // printf("sum2= %.10lf\n",sum/(DL*DMenergy(m,r,dr)));
  printf("execution time: %lf\n",time_spent);

  // double low=fabs(m[2]);
  // int j=0;
  for (i=0;i<3*N;i++){
      printf("%lf ",m[i]);

    }

  // printf("j: %d \n",j);
  // printf("m_z low = %lf \n",low);


  free(t);
  free(y);
  return 0;
}



float second_derivative(double *m,int i,double dr){
  return (m[i+3]+m[i-3]-2*m[i])/pow(dr,2);

}

float first_derivative(double *m,int i,double dr){
  return (m[i+3]-m[i-3])/(2*dr);
}

void laplace(double *m,double *r,double *laplacian,int i,double dr){
  int j=i/3;
  laplacian[0]=second_derivative(m,i,dr)+first_derivative(m,i,dr)/r[j]-m[i]/pow(r[j],2);
  laplacian[1]=second_derivative(m,i+1,dr)+first_derivative(m,i+1,dr)/r[j]-m[i+1]/pow(r[j],2);
  laplacian[2]=second_derivative(m,i+2,dr)+first_derivative(m,i+2,dr)/r[j];
}

void skyrmion_init_values(double *m0,double *r){
  int i,j;

  m0[0]=0;m0[1]=0;m0[2]=-1;
  m0[3*N-3]=0;m0[3*N-2]=0;m0[3*N-1]=1;

  for (i=3;i<3*N-3;i=i+3){
    j=i/3;
    // ---- Domain Wall -----
    //
    // m0[i]=1/cosh(r[j]-5);
    // m0[i+2]=tanh(r[j]-5);

    //---- Belavin Polyakov ----
    m0[i+1]=2*r[j]/(1+r[j]*r[j]);
    m0[i+2]=(r[j]*r[j]-1)/(r[j]*r[j]+1);

    //---- Skyrmionium ----

    //--- transform----
    // double n1=1/cosh(r[j]-3);
    // double n3 =tanh(r[j]-3);
    //
    // m0[i+0]=2*n3*n1;
    // m0[i+2]=2*n3*n3-1;

    // ---ansatz---
    // double C=(r[j]*r[j]+1)/r[j];
    // m0[i+1]=2*C/(1+C*C);
    // m0[i+2]=(1-C*C)/(1+C*C);
  }
}

double ANenergy(double *m, double *r, double dr){
  double integral = 0;
  for (int i=0;i<3*N;i=i+3){
    int j = i/3;
    integral+=(1-pow(m[i+2],2))*r[j];
  }
  integral*=M_PI*dr;
  return integral;
}

double EXenergy(double *m, double *r, double dr){
  double intergal = 0;
  for (int i=0;i<3*N;i=i+3){
    int j=i/3;
    if (r[j]==r[0] || r[j]==r[N-1]){
      continue;
    }
    intergal+=(pow(first_derivative(m,i,dr),2)+pow(first_derivative(m,i+1, dr),2)\
            +pow(first_derivative(m,i+2, dr),2)+(pow(m[i],2)+pow(m[i+1],2))/r[j])*r[j];


  }
  intergal*=M_PI*dr;
  return intergal;
}

double DMenergy(double *m, double *r, double dr){
  double integral = 0;
  for (int i=0;i<3*N;i=i+3){
    int j = i/3;
    if (r[j]==r[0] || r[j]==r[N-1]){
      continue;
    }
    integral+=(first_derivative(m, i+1, dr)*m[i+2]-first_derivative(m, i+2, dr)*m[i+1]+m[i+1]*m[i+2]/r[j])*r[j];
  }
  integral*= 2*DL*M_PI*dr;
  return integral;
}




void LLequation(double t, double *m,double *dm,double *r,double dr){
  float crProd[3]={0};
  float dump[3]={0};
  double f[3]={0};
  int i,j,k;
  double dotProd;

  // for(i=0;i<3*N;i++){
  //   printf("m %f\n",m[i]);
  // }

  for(i=3;i<3*N-3;i=i+3){
    dotProd=0;
    j=i/3;

    laplace(m,r,f,i,dr);
    f[1]+=2*DL*(first_derivative(m,i+2,dr));
    f[2]+=-2*DL*(first_derivative(m,i+1,dr)+m[i+1]/r[j]);
    f[2]+=m[i+2];
    if(f[0]!=f[0] || f[1]!=f[1] || f[2]!=f[2]){
      printf("problem with f\n");
      exit(1);
    }
    // cross product
    crProd[0]=m[i+1]*f[2]-m[i+2]*f[1];
    crProd[1]=m[i+2]*f[0]-m[i]*f[2];
    crProd[2]=m[i]*f[1]-m[i+1]*f[0];
    if(crProd[0]!=crProd[0] || crProd[1]!=crProd[1] || crProd[2]!=crProd[2]){
      printf("problem with crProd\n");
      exit(1);
    }    // dot product
    for(k=0;k<3;k++){
      dotProd+=m[i+k]*f[k];
    }
    dump[0]=DP*(dotProd*m[i]-f[0]);
    dump[1]=DP*(dotProd*m[i+1]-f[1]);
    dump[2]=DP*(dotProd*m[i+2]-f[2]);
    if(dump[0]!=dump[0] || dump[1]!=dump[1] || dump[2]!=dump[2]){
      printf("problem with dump\n");
      exit(1);
    }
    dm[i]=-crProd[0]-dump[0];
    dm[i+1]=-crProd[1]-dump[1];
    dm[i+2]=-crProd[2]-dump[2];
    // if(dm[i]!=dm[i] || dm[i+1]!=dm[i+1] || dm[i+2]!=dm[i+2]){
    //   print()
    // }

  }


}
