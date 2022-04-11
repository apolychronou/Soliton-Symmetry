#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <complex.h>
#include "rk4.h"
#include "space.h"
#define N 400 //r grid
#define NZ 200 //z grid
#define DL 0.63 //d-m interaction coef
// #define DL 0.0
#define DP 0.005//dumping coef
#define AN 1  //anisotropy coef
#define APROJ 1 //stereographic projection
#define PHI M_PI // stereographic projection angle
#define T_STOP 0
#define INIT 1


float deriv2r(double *m,int i,double dr);
float deriv1r(double *m,int i, double dr);
float deriv2z(double *m,int i,double dz);
float deriv1r(double *m,int i,double dz);
void laplace(double *m,double *r,double *z,double *laplacian,int ri,int zi,double dr,double dz);
void vortexRing_init_values(double *m0,double *r, double *z,double a,double phi);
void LLequation(double t, double *m,double *dm,double *r,double *z,double dr,double dz);
void boundary_conditions(double *m);
double EXenergy(double *m, double *r,double *z, double dr,double dz);
double ANenergy(double *m, double *r,double *z, double dr,double dz);
double DMenergy(double *m, double *r,double *z, double dr,double dz);
void DMenergy_Density(double *m, double *r,double *z, double dr,double dz);



int main(int argc, char *argv[]){
  // double *y;
  // double m0[3*N*NZ]={0};
  double *m0;
  int i=0;
  double dr=0.1;
  double dz=0.1;
  double r[N];
  double z[NZ];
  int steps;
  double *t;
  double tspan[2];
  double a,phi=0;
  FILE *f;
  // double m[3*N*NZ]={0};
  // double prevm[3*N*NZ];
  // double y[3*N*NZ*2]={0};
  double *y,*m,*prevm;
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
  m=(double*)malloc((3*N*NZ)*sizeof(double));
  prevm=(double*)malloc((3*N*NZ)*sizeof(double));
  m0=y;

  vortexRing_init_values(y,r,z,1,phi);

  DMenergy_Density(y,r,z,dr,dz);
  return(0);

  if(INIT){
    f=fopen("./data/center.csv","w+");
    fprintf(f,"a,phi,E_ex,E_an,E_dm,Sum\n");

    for(a=dr;a<dr*N;a+=dr){
      vortexRing_init_values(y,r,z,a,phi);
      double init_dm=DMenergy(y,r,z,dr,dz);
      double init_ex=EXenergy(y,r,z,dr,dz);
      double init_an=ANenergy(y,r,z,dr,dz);
      fprintf(f,"%lf, %lf, %lf, %lf, %lf, %lf\n",a,phi,init_ex,
      init_an,init_dm,init_ex+2*init_dm+3*init_an);
    }
    fclose(f);

    f=fopen("./data/angle.csv","w+");
    fprintf(f,"a,phi,E_ex,E_an,E_dm,Sum\n");
    a=1;
    for(phi=0;phi<M_PI;phi+=0.01){
      vortexRing_init_values(y,r,z,a,phi);
      double init_dm=DMenergy(y,r,z,dr,dz);
      double init_ex=EXenergy(y,r,z,dr,dz);
      double init_an=ANenergy(y,r,z,dr,dz);
      fprintf(f,"%lf, %lf, %lf, %lf, %lf, %lf\n",a,phi,init_ex,
      init_an,init_dm,init_ex+2*init_dm+3*init_an);
    }

    // for(i=0;i<3*N*NZ;i++){
    //   printf("%.5lf ",m0[i]);
    // }
    return 1;
  }


  clock_t begin = clock();
  rk4 ( LLequation, tspan, m0, steps, 3*N*NZ, t, y,r,z,dr,dz );

  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  //
  // double ANener = ANenergy(m,r,z,dr,dz);
  // double DMener = DMenergy(m,r,z,dr,dz);
  // double EXener = EXenergy(m,r,z,dr,dz);
  // double sum=2*DMener+3*ANener+EXener;
  // printf("Energy: %.10lf\t", ANener+DMener+EXener);
//   printf("Prev energy: %.10lf\n",ANenergy(prevm,r,z,dr,dz)+ DMenergy(prevm,r,z,dr,dz) + EXenergy(prevm,r,z,dr,dz));
  // printf("Initial energy: %.10lf\t",ANenergy(m0,r,z,dr,dz)+ DMenergy(m0,r,z,dr,dz) + EXenergy(m0,r,z,dr,dz));
//   if(ANener!=0){
//   printf("sum1= %.10lf\n",sum/EXener);
//   }
  printf("execution time: %lf\n",time_spent);

  // for (i=0;i<3*N*NZ;i++){
  //   // m[i]=y[1*3*N*NZ+i];
  //   printf("%.5lf ",y[1*3*N*NZ+i]);
  //
  // }
//   int k=0;
//   for (i=0;i<3*N*NZ;i=i+3*N){
//   for (k=0;k<3*N;k++){
//       printf("%lf ",m[i+k]);
//     }
//   }


  free(t);
  free(y);
  free(m);
  free(prevm);
  return 0;
}



float deriv2r(double *m,int i,double dr){
  return (m[i+3]+m[i-3]-2*m[i])/(dr*dr);

}

float deriv1r(double *m,int i,double dr){
  return (m[i+3]-m[i-3])/(2*dr);
}

float deriv2z(double *m,int i,double dz){
  return (m[i+3*N]+m[i-3*N]-2*m[i])/(dz*dz);
}

float deriv1z(double *m,int i ,double dz){
  return (m[i+3*N]-m[i-3*N])/(2*dz);
}

void laplace(double *m,double *r,double *z,double *laplacian,int ri, int zi,double dr,double dz){
  int j=ri/3;

  laplacian[0]=deriv2r(m,zi+ri,dr)+deriv1r(m,zi+ri,dr)/r[j]-m[ri+zi]/(r[j]*r[j])+
               deriv2z(m,zi+ri,dz);
  laplacian[1]=deriv2r(m,ri+zi+1,dr)+deriv1r(m,ri+zi+1,dr)/r[j]+
               deriv2z(m,zi+ri+1,dz)-m[ri+zi+1]/(r[j]*r[j]);
  laplacian[2]=deriv2r(m,ri+zi+2,dr)+deriv1r(m,ri+zi+2,dr)/r[j]+
               +deriv2z(m,zi+ri+2,dz);
}

void vortexRing_init_values(double *m0,double *r,double *z,double a,double phi){
  int ri,zi,i,j;
  // double complex omega,omega_bar;
  // double a=1*APROJ;
  // double A,B,C,E;
  double real,im;

  for (zi=0;zi<3*N*NZ;zi=zi+3*N){
    for (ri=0;ri<3*N;ri=ri+3){
      i=ri/3;
      j=zi/3/N;
      if (r[i]==r[0] || r[i]==r[N-1] || z[j]==z[NZ-1] || z[j]==z[0]){
        m0[ri+zi]=0;m0[ri+zi+1]=0;m0[ri+zi+2]=1;
      }
      else if(fabs(r[i]-1*a)<1e-10 && fabs(z[j])<1e-10){
        m0[ri+zi]=0;m0[ri+zi+1]=0;m0[ri+zi+2]=-1;
      }
      else{
        /*
        // double c=2*a*r[i]/(pow(2*a*z[j],2)+pow(r[i]*r[i]+z[j]*z[j]-a*a,2));
        // m0[ri+zi+0]=-c*4*a*z[j]/(1+c*2*a*r[i]);
        // m0[ri+zi+1]=c*2*(r[i]*r[i]+z[j]*z[j]-a*a)/(1+c*2*a*r[i]);
        // m0[ri+zi+2]=(1-c*2*a*r[i])/(1+c*2*a*r[i]);

        // double scalen=exp(-(r[i]+z[j])/a);
        // double scalep=exp((r[i]+z[j])/a);
        // scalen=1;scalep=1;
        // A=2*a*z[j];B=r[i]*r[i]+z[j]*z[j]-a*a;
        // C=(A*A+B*B)*scalep+scalen*(2*a*r[i])*(2*a*r[i]);
        // m0[ri+zi+0]=2*a*r[i]*(2*A*cos(PHI*1)+2*B*sin(PHI*1))/C;
        // m0[ri+zi+1]=2*a*r[i]*(2*A*sin(PHI*1)-2*B*cos(PHI*1))/C;
        // m0[ri+zi+2]=((A*A+B*B)*scalep-(2*a*r[i])*(2*a*r[i])*scalen)/C;

        // A=2*a*z[j];B=r[i]*r[i]+z[j]*z[j]-a*a;
        // C=2*a*r[i]/(pow(2*a*z[j],2)+pow(r[i]*r[i]+z[j]*z[j]-a*a,2));
        // // E=exp(-r[i]*r[i]/a*a)*exp(-z[j]*z[j]/a*a);
        // E=1;
        // C=C*E;
        // m0[ri+zi+0]=(2*A*cos(PHI*1)+2*B*sin(PHI*1))*C/(1+C*C*pow(r[i]*r[i]+z[j]*z[j]-a*a,2));
        // m0[ri+zi+1]=(2*A*sin(PHI*1)-2*B*cos(PHI*1))*C/(1+C*C*pow(r[i]*r[i]+z[j]*z[j]-a*a,2));
        // m0[ri+zi+2]=(1-C*C*pow(r[i]*r[i]+z[j]*z[j]-a*a,2))/(1+C*C*pow(r[i]*r[i]+z[j]*z[j]-a*a,2));



        // A=2*a*z[j];B=r[i]*r[i]+z[j]*z[j]-a*a;
        // C=pow(2*a*r[i],2)/(A*A+B*B+pow(2*a*r[i],2));
        // m0[ri+zi+0]=C*(2*A*cos(PHI*1)+2*B*sin(PHI*1));
        // m0[ri+zi+1]=C*(2*A*sin(PHI*1)-2*B*cos(PHI*1));
        // m0[ri+zi+2]=(A*A+B*B-pow(2*a*r[i],2))/(A*A+B*B+pow(2*a*r[i],2));

        */
        double C,E;
        C=2*a*r[i]/(pow(2*a*z[j],2)+pow(r[i]*r[i]+z[j]*z[j]-a*a,2));
        E=exp(-r[i]*r[i]/(a*a))*exp(-(z[j]*z[j])/(a*a));
        // E=exp(-r[i]/(a))*exp(-fabs(z[j])/(a));
        E=1;
        C=C*E;
        // real=+2*a*z[j]*C;
        // im=-(r[i]*r[i]+z[j]*z[j]-a*a)*C;
        real=(2*a*z[j]*cos(phi)+(r[i]*r[i]+z[j]*z[j]-a*a)*sin(phi))*C;
        im=(2*a*z[j]*sin(phi)-(r[i]*r[i]+z[j]*z[j]-a*a)*cos(phi))*C;

        double oobar=(real*real+im*im);

        m0[ri+zi+0]=(2*real)/(1+oobar);
        m0[ri+zi+1]=(2*im)/(1+oobar);
        m0[ri+zi+2]=(1-oobar)/(1+oobar);
      }

    }
  }

}

void skyrmion_init_values(double *m0,double *r,double *z){
  int i,ri,zi,j;
  // int k=5;
  //
  // k=2;
  for (zi=0;zi<3*N*NZ;zi=zi+3*N){
    for (ri=0;ri<3*N;ri=ri+3){
        i=ri/3;
        j=zi/3/N;
        if (r[i]==r[0]){
          m0[ri+zi]=0;m0[ri+zi+1]=0;m0[ri+zi+2]=-1;
        }else if (r[i]==r[N-1]){
          m0[ri+zi]=0;m0[ri+zi+1]=0;m0[ri+zi+2]=1;
        // }else if (z[j]==z[0]||z[j]==z[NZ-1]) {
        //   m0[ri+zi]=0;m0[ri+zi+1]=0;m0[ri+zi+2]=1;
      }else if (z[j]<0){
          // ---------DOMAIN WALL------
          // m0[ri+zi]=1/cosh(r[i]-2);
          // m0[ri+zi+2]=tanh(r[i]-2);

          //---------BELAVIN POLYAKOV----
          m0[ri+zi]=2*(r[i])/(1+(r[i])*(r[i]));
          m0[ri+zi+2]=((r[i])*(r[i])-1)/((r[i])*(r[i])+1);
        }
        else{
          // ---------DOMAIN WALL------
          // m0[ri+zi]=1/cosh(r[i]-6);
          // m0[ri+zi+2]=tanh(r[i]-6);

          //---------BELAVIN POLYAKOV----
          m0[ri+zi+1]=2*r[i]/(1+r[i]*r[i]);
          m0[ri+zi+2]=(r[i]*r[i]-1)/(r[i]*r[i]+1);
        }
      }
      // k=((k+1)%15!=0)?k+1:5;
      // k=((k%4)==0)?2:k+1;
    }
}

double ANenergy(double *m, double *r,double *z, double dr,double dz){
  double integral = 0;
  int zi,ri,i;
  for (zi=0;zi<3*N*NZ;zi=zi+3*N){
    for (ri=0;ri<3*N;ri=ri+3){
       i = ri/3;
      integral+=(1-pow(m[ri+zi+2],2))*r[i];
    }
  }
  integral*=M_PI*dr*dz;
  return integral;
}

double EXenergy(double *m, double *r,double *z, double dr,double dz){
  double integral = 0;
  int i,j,ri,zi;

  for (zi=0;zi<3*N*NZ;zi=zi+3*N){
    for(ri=0;ri<3*N;ri=ri+3){
      i=ri/3;j=zi/3/N;
      if (r[i]==r[0] || r[i]==r[N-1] || z[j]==z[0] || z[j]==z[NZ-1]){
        continue;
      }
        integral+=(pow(deriv1r(m,ri+zi,dr),2)+pow(deriv1r(m,zi+ri+1, dr),2)
                +pow(deriv1r(m,zi+ri+2, dr),2)+(pow(m[zi+ri],2)+pow(m[zi+ri+1],2))/r[i]+
                pow(deriv1z(m,zi+ri,dz),2)+pow(deriv1z(m,zi+ri+1, dz),2)+
                pow(deriv1z(m,zi+ri+2, dz),2))*r[i];
        if (isinf(integral)){
          printf("zi:%d ri:%d i:%d j:%d",zi,ri,i,j);
          break;
        }


    }
  }
  integral*=M_PI*dr*dz;
  return integral;
}

double DMenergy(double *m, double *r,double *z, double dr,double dz){
  double integral = 0;
  int zi,ri,i;
  for (zi=0;zi<3*N*NZ;zi=zi+3*N){
    for (ri=0;ri<3*N;ri=ri+3){
      i = ri/3;
      // j = zi/3/N;

      // if (r[i]==r[0] || r[i]==r[N-1] || z[j]==z[0] || z[j]==z[NZ-1]){
      if (r[i]==r[0] || r[i]==r[N-1] ){
        continue;
        }
        else{
          integral+=(deriv1r(m, ri+zi+1, dr)*m[ri+zi+2]-deriv1r(m, zi+ri+2, dr)*m[zi+ri+1]+m[zi+ri+1]*m[zi+ri+2]/r[i])*r[i];
        }
      }
  }
  integral*= 2*DL*M_PI*dr*dz;
  return integral;
}

void boundary_conditions(double *m){
  int i;
  for(i=0;i<3*N;i++){
    m[i]=m[i+3*N];
    m[3*N*(NZ-1)+i]=m[3*N*(NZ-2)+i];
  }
  return;
}


void LLequation(double t, double *m,double *dm,double *r,double *z,double dr, double dz){
  float crProd[3]={0};
  float dump[3]={0};
  double f[3]={0};
  int ri,zi,j,k;
  double dotProd;

  // for(i=0;i<3*N;i++){
  //   printf("m %f\n",m[i]);
  // }
  for (zi=3*N;zi<3*N*(NZ-1);zi=zi+3*N){
    for(ri=3;ri<3*N-3;ri=ri+3){
      dotProd=0;
      j=ri/3;


      laplace(m,r,z,f,ri,zi,dr,dz);
      f[1]+=2*DL*(deriv1r(m,ri+zi+2,dr));
      f[2]+=-2*DL*(deriv1r(m,ri+zi+1,dr)+m[ri+zi+1]/r[j]);
      f[2]+=AN*m[ri+zi+2];
      if(f[0]!=f[0] || f[1]!=f[1] || f[2]!=f[2]){
        printf("problem with f\n");
        exit(1);
      }
      // cross product
      crProd[0]=m[ri+zi+1]*f[2]-m[ri+zi+2]*f[1];
      crProd[1]=m[ri+zi+2]*f[0]-m[ri+zi]*f[2];
      crProd[2]=m[ri+zi]*f[1]-m[ri+zi+1]*f[0];
      if(crProd[0]!=crProd[0] || crProd[1]!=crProd[1] || crProd[2]!=crProd[2]){
        printf("problem with crProd\n");
        exit(1);
      }
      // dot product
      for(k=0;k<3;k++){
        dotProd+=m[ri+zi+k]*f[k];
      }
      dump[0]=DP*(dotProd*m[ri+zi]-f[0]);
      dump[1]=DP*(dotProd*m[ri+zi+1]-f[1]);
      dump[2]=DP*(dotProd*m[ri+zi+2]-f[2]);
      if(dump[0]!=dump[0] || dump[1]!=dump[1] || dump[2]!=dump[2]){
        printf("problem with dump\n");
        exit(1);
      }
      dm[ri+zi]=-crProd[0]-dump[0];
      dm[ri+zi+1]=-crProd[1]-dump[1];
      dm[ri+zi+2]=-crProd[2]-dump[2];

    }
  }
  // boundary_conditions(dm);
  // ;
}


void DMenergy_Density(double *m, double *r,double *z, double dr,double dz){
  double integral =0;
  int zi,ri,i,j;
  // for (zi=NZ/2;zi<NZ/2+1;zi++){
    printf("z,r,E_dm\n");
    for (ri=0;ri<3*N;ri=ri+3){
      i = ri/3;
      j=NZ/2;
      zi = j*3*N;

      // if (r[i]==r[0] || r[i]==r[N-1] || z[j]==z[0] || z[j]==z[NZ-1]){
      if (r[i]==r[0] || r[i]==r[N-1] ){
        continue;
        }
        else{
          integral=(deriv1r(m, ri+zi+1, dr)*m[ri+zi+2]-deriv1r(m, zi+ri+2, dr)*m[zi+ri+1]+m[zi+ri+1]*m[zi+ri+2]/r[i])*r[i];
          printf("%lf,%lf,%lf\n",z[j],r[i],integral);
        }
      }
  // }
  integral*= 2*DL*M_PI*dr*dz;
  // return integral;
}

/*
void DMenergy_Density(double *m, double *r,double *z, double dr,double dz){
  double integral =0;
  int zi,ri,i,j;
  for (zi=NZ/2;zi<NZ/2+1;zi++){
    printf("z,r,E_dm\n");
    for (ri=0;ri<3*N;ri=ri+3){
      i = ri/3;
      // j=NZ/2;
      // zi = j*3*N;
      j=zi/3/N;

      // if (r[i]==r[0] || r[i]==r[N-1] || z[j]==z[0] || z[j]==z[NZ-1]){
      if (r[i]==r[0] || r[i]==r[N-1] ){
        continue;
        }
        else{
          integral=(deriv1r(m, ri+zi+1, dr)*m[ri+zi+2]-deriv1r(m, zi+ri+2, dr)*m[zi+ri+1]+m[zi+ri+1]*m[zi+ri+2]/r[i])*r[i];
          printf("%lf,%lf,%lf\n",z[j],r[i],integral);
        }
      }
  }
  integral*= 2*DL*M_PI*dr*dz;
  // return integral;
}
*/
