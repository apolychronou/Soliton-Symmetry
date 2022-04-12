#include "axialSym.h"
// #include <stdio.h>
// #include <string.h>
// #include <math.h>
// #include <stdlib.h>
// #define N 400 //r grid
// #define NZ 200 //z grid
// #define DL 0.63 //d-m interaction coef
// #define DP 0.005//dumping coef
// #define AN 1  //anisotropy coef

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

void laplace(double *m,double *r,double *z,double *laplacian,int i, int j,double dr,double dz){
  int ind=magnIndex(i,j);
  laplacian[0]=deriv2r(m,ind,dr)+deriv1r(m,ind,dr)/r[i]-m[ind]/(r[i]*r[i])+
               deriv2z(m,ind,dz);
  laplacian[1]=deriv2r(m,ind+1,dr)+deriv1r(m,ind+1,dr)/r[i]+
               deriv2z(m,ind+1,dz)-m[ind+1]/(r[i]*r[i]);
  laplacian[2]=deriv2r(m,ind+2,dr)+deriv1r(m,ind+2,dr)/r[i]+
               +deriv2z(m,ind+2,dz);
}

void vortexRing_init_values(double *m0,double *r,double *z,double a,double phi){
  int i=0,j=0,ind=0;
  double real=0,im=0;

  for (j=0;j<NZ;j++){
    for (i=0;i<N;i++){
      ind=magnIndex(i,j);
      if (i==0 || i==N-1 || j==NZ-1 || j==0){
        m0[ind]=0;m0[ind+1]=0;m0[ind+2]=1;
      }
      else if(fabs(r[i]-a)<1e-10 && fabs(z[j])<1e-10){
        m0[ind]=0;m0[ind+1]=0;m0[ind+2]=-1;
      }
      else{
        double C,E;
        C=2*a*r[i]/(pow(2*a*z[j],2)+pow(r[i]*r[i]+z[j]*z[j]-a*a,2));
        E=exp(-r[i]*r[i]/(a*a))*exp(-(z[j]*z[j])/(a*a));
        // E=exp(-r[i]/(a))*exp(-fabs(z[j])/(a));
        // E=1;
        C=C*E;
        real=(2*a*z[j]*cos(phi)+(r[i]*r[i]+z[j]*z[j]-a*a)*sin(phi))*C;
        im=(2*a*z[j]*sin(phi)-(r[i]*r[i]+z[j]*z[j]-a*a)*cos(phi))*C;

        double oobar=(real*real+im*im);

        m0[ind+0]=(2*real)/(1+oobar);
        m0[ind+1]=(2*im)/(1+oobar);
        m0[ind+2]=(1-oobar)/(1+oobar);
      }

    }
  }

}

int magnIndex(int i, int j){
  int ind=0;
  ind=3*N*j+3*i;
  return ind;
}

void LLequation(double t, double *m,double *dm,double *r,double *z,double dr, double dz){
  double crProd[3]={0};
  double damp[3]={0};
  double f[3]={0};
  int j=0,i=0,k=0;
  double dotProd=0;
  int ind=0;

  for (j=1;j<NZ-1;j++){
      for(i=1;i<N-1;i++){
      dotProd=0;
      ind=magnIndex(i,j);

      laplace(m,r,z,f,i,j,dr,dz);
      f[1]+=2*DL*(deriv1r(m,ind+2,dr));
      f[2]+=-2*DL*(deriv1r(m,ind+1,dr)+m[ind+1]/r[i]);
      f[2]+=AN*m[ind+2];
      if(f[0]!=f[0] || f[1]!=f[1] || f[2]!=f[2]){
        printf("problem with f\n");
        exit(1);
      }
      // cross product
      crProd[0]=m[ind+1]*f[2]-m[ind+2]*f[1];
      crProd[1]=m[ind+2]*f[0]-m[ind]*f[2];
      crProd[2]=m[ind]*f[1]-m[ind+1]*f[0];
      if(crProd[0]!=crProd[0] || crProd[1]!=crProd[1] || crProd[2]!=crProd[2]){
        printf("problem with crProd\n");
        exit(1);
      }
      // dot product
      for(k=0;k<3;k++){
        dotProd+=m[ind+k]*f[k];
      }
      damp[0]=DP*(dotProd*m[ind]-f[0]);
      damp[1]=DP*(dotProd*m[ind+1]-f[1]);
      damp[2]=DP*(dotProd*m[ind+2]-f[2]);
      if(damp[0]!=damp[0] || damp[1]!=damp[1] || damp[2]!=damp[2]){
        printf("problem with damp\n");
        exit(1);
      }
      dm[ind]=-crProd[0]-damp[0];
      dm[ind+1]=-crProd[1]-damp[1];
      dm[ind+2]=-crProd[2]-damp[2];

    }
  }
  // boundary_conditions(dm);
  // ;
}

double DMenergy(double *m, double *r,double *z, double dr,double dz){
  double integral = 0;
  int ind,i,j;
  for(j=0;j<NZ;j++){
    for (i=0;i<N;i++){
      ind=magnIndex(i,j);

      if (i==0 || i==N-1 ){
        continue;
        }
        else{
          integral+=(deriv1r(m, ind+1, dr)*m[ind+2]-deriv1r(m, ind+2, dr)*m[ind+1]+m[ind+1]*m[ind+2]/r[i])*r[i];
        }
      }
  }
  integral*= 2*DL*M_PI*dr*dz;
  return integral;
}

double EXenergy(double *m, double *r,double *z, double dr,double dz){
  double integral = 0;
  int i=0,j=0,ind=0;

  for (j=0;j<NZ;j++){
    for(i=0;i<N;i++){
      ind=magnIndex(i,j);
      if (i==0 || i==N-1 || j==0 || j==NZ-1){
        continue;
      }
        integral+=(pow(deriv1r(m,ind,dr),2)+pow(deriv1r(m,ind+1, dr),2)
                +pow(deriv1r(m,ind+2, dr),2)+(pow(m[ind],2)+pow(m[ind+1],2))/r[i]+
                pow(deriv1z(m,ind,dz),2)+pow(deriv1z(m,ind+1, dz),2)+
                pow(deriv1z(m,ind+2, dz),2))*r[i];
    }
  }
  integral*=M_PI*dr*dz;
  return integral;
}

double ANenergy(double *m, double *r,double *z, double dr,double dz){
  double integral = 0;
  int j=0,i=0,ind=0;
  for (j=0;j<NZ;j++){
    for (i=0;i<N;i++){
      ind=magnIndex(i,j);
      integral+=(1-pow(m[ind+2],2))*r[i];
    }
  }
  integral*=M_PI*dr*dz;
  return integral;
}

void DMenergy_Density(double *m, double *r,double *z, double dr,double dz,double a,double phi){
  double integral =0;
  int ind,i,j;
  FILE *f;
  vortexRing_init_values(m,r,z,a,phi);
  f=fopen("./data/density.csv","w+");

    j=NZ/2;
    fprintf(f,"z,r,E_dm\n");
    for (i=0;i<N;i++){
      ind=magnIndex(i,j);

      if (i==0 || i==N-1 ){
        continue;
        }
        else{
          integral=(deriv1r(m, ind+1, dr)*m[ind+2]-deriv1r(m, ind+2, dr)*m[ind+1]+m[ind+1]*m[ind+2]/r[i])*r[i];
          fprintf(f,"%lf,%lf,%lf\n",z[j],r[i],integral);
        }
      }
    fclose(f);
}

void DMenergy_zDensity(double *m, double *r,double *z, double dr,double dz,double a,double phi){
  double integral =0;
  int ind=0,i=0,j=0;
  FILE *f;
  vortexRing_init_values(m,r,z,a,phi);
  f=fopen("./data/z_density.csv","w+");
    fprintf(f,"z,r,E_dm\n");
    for(j=NZ/2;j<NZ;j++){
      integral=0;
      for (i=0;i<N;i++){
        ind=magnIndex(i,j);
        if (i==0 || i==N-1 ){
          continue;
        }
        else{
            integral+=(deriv1r(m, ind+1, dr)*m[ind+2]-deriv1r(m, ind+2, dr)*m[ind+1]+m[ind+1]*m[ind+2]/r[i])*r[i];
        }
      }
      fprintf(f,"%lf,%lf,%lf\n",z[j],r[i],integral);
    }
    fclose(f);
}

void center_angle(double *y,double *r,double *z, double dr, double dz){
  double a=0,phi=0;
  FILE *f=NULL;
  phi=0;

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

}

void skyrmion_init_values(double *m0,double *r,double *z){
  int i=0,ind=0,j=0;
  // double n1=0,n3=0;
  //
  // k=2;
  for (j=0;j<NZ;j++){
    for (i=0;i<N;i++){
      ind=magnIndex(i,j);
        if (i==0){
          m0[ind]=0;m0[ind+1]=0;m0[ind+2]=-1;
        }else if (i==N-1){
          m0[ind]=0;m0[ind+1]=0;m0[ind+2]=1;
        // }else if (z[j]==z[0]||z[j]==z[NZ-1]) {
        //   m0[ind]=0;m0[ind+1]=0;m0[ind+2]=1;
      }else if (z[j]<0){
          // ---------DOMAIN WALL------
          // m0[ind]=1/cosh(r[i]-2);
          // m0[ind+2]=tanh(r[i]-2);

          //-----------SKRYRMIONIUM TRANSFORMATION
          // n1=1/cosh(r[i]-3);
          // n3=tanh(r[i]-3);
          // m0[ind+1]=2*n3*n1;
          // m0[ind+2]=2*n3*n3-1;


          //---------BELAVIN POLYAKOV----
          m0[ind]=2*(r[i])/(1+(r[i])*(r[i]));
          m0[ind+1]=0;
          m0[ind+2]=((r[i])*(r[i])-1)/((r[i])*(r[i])+1);

        }else{
          // ---------DOMAIN WALL------
          // m0[ind]=1/cosh(r[i]-6);
          // m0[ind+2]=tanh(r[i]-6);


          //-----------SKRYRMIONIUM TRANSFORMATION
          // n1=1/cosh(r[i]-2);
          // n3 =tanh(r[i]-2);
          // m0[ind+1]=2*n3*n1;
          // m0[ind+2]=2*n3*n3-1;


          //---------BELAVIN POLYAKOV----
          m0[ind+1]=2*r[i]/(1+r[i]*r[i]);
          m0[ind+0]=0;
          m0[ind+2]=(r[i]*r[i]-1)/(r[i]*r[i]+1);
        }
      }

    }
}

void boundary_conditions(double *m){
  int i=0;
  for(i=3;i<3*N-3;i++){
    m[i]=m[i+3*N];
    m[3*N*(NZ-1)+i]=m[3*N*(NZ-2)+i];
  }
  return;
}

void moving_LLequation(double t, double *m,double *dm,double *r,double *z,double dr, double dz){
  double crProd[3]={0};
  double damp[3]={0};
  double f[3]={0};
  double deriv[3]={0};
  int j=0,i=0,k=0;
  double dotProd=0;
  int ind=0;
  double v=V;



  for (j=1;j<NZ-1;j++){
      for(i=1;i<N-1;i++){
      dotProd=0;
      ind=magnIndex(i,j);

      laplace(m,r,z,f,i,j,dr,dz);
      f[1]+=2*DL*(deriv1r(m,ind+2,dr));
      f[2]+=-2*DL*(deriv1r(m,ind+1,dr)+m[ind+1]/r[i]);
      f[2]+=AN*m[ind+2];
      if(f[0]!=f[0] || f[1]!=f[1] || f[2]!=f[2]){
        printf("problem with f\n");
        exit(1);
      }
      deriv[0]=(deriv1z(m,ind,dz))*v;
      deriv[1]=(deriv1z(m,ind+1,dz))*v;
      deriv[2]=(deriv1z(m,ind+2,dz))*v;


      // cross product
      crProd[0]=m[ind+1]*deriv[2]-m[ind+2]*deriv[1];
      crProd[1]=m[ind+2]*deriv[0]-m[ind]*deriv[2];
      crProd[2]=m[ind]*deriv[1]-m[ind+1]*deriv[0];
      if(crProd[0]!=crProd[0] || crProd[1]!=crProd[1] || crProd[2]!=crProd[2]){
        printf("problem with crProd\n");
        exit(1);
      }
      // dot product
        dotProd+=m[ind+0]*f[0];
        dotProd+=m[ind+1]*f[1];
        dotProd+=m[ind+2]*f[2];

      damp[0]=(dotProd*m[ind]-f[0]);
      damp[1]=(dotProd*m[ind+1]-f[1]);
      damp[2]=(dotProd*m[ind+2]-f[2]);
      if(damp[0]!=damp[0] || damp[1]!=damp[1] || damp[2]!=damp[2]){
        printf("problem with damp\n");
        exit(1);
      }
      dm[ind]=crProd[0]-damp[0];
      dm[ind+1]=crProd[1]-damp[1];
      dm[ind+2]=crProd[2]-damp[2];

    }
  }
  // boundary_conditions(dm);
  // ;
}

double mz_integral(double *m,double *r,double *z, double dr, double dz){
  double integral = 0;
  int j,i,ind;
  for (j=0;j<NZ;j++){
    for (i=0;i<N;i++){
      ind=magnIndex(i,j);
      integral+=(1-m[ind+2])*r[i];
    }
  }
  integral*=M_PI*dr*dz;
  return integral;
}
