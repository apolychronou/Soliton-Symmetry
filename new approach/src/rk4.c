// # include <math.h>
// # include <stdio.h>
// # include <stdlib.h>
// # include <time.h>
// #include <string.h>
// #include <errno.h>
# include "axialSym.h"
// # include "rk4.h"

/******************************************************************************/

void rk4 ( void dydt ( double t, double u[], double f[], double *r,double *z,double dr,double dz ), double tspan[2],
  double y0[], int n, int m, double t[], double y[], double *r,double *z, double dr,double dz )

/******************************************************************************/
/*
  Purpose:

    rk4 approximates an ODE using a Runge-Kutta fourth order method.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 April 2020

  Author:

    John Burkardt

  Input:

    double DYDT ( double T, double U ), a function which evaluates
    the derivative, or right hand side of the problem.

    double TSPAN[2]: the initial and final times

    double Y0[M]: the initial condition

    int N: the number of steps to take.

    int M: the number of variables.

  Output:

    double t[n+1], y[(n+1)*m]: the times and solution values.


  Disclaimer: The above code was the basis and is changed in order to suite
  this project's needs
*/
{
  double dt;
  double *f0;
  double *f1;
  double *f2;
  double *f3;
  int i;
  int j;
  double t0;
  double t1;
  double t2;
  double t3;
  double *u0;
  double *u1;
  double *u2;
  double *u3;
  double e_ex=0,e_dm=0,e_an=0,time_point;
  double en=0;double preven=3e20;
  double norm=0;
  int points,tpoints;
  int filecounter=0;
  char *countername;
  char filename[42];
  FILE *f=NULL;
  extern int errno;
  errno=0;
  int errnum;
  int prfiles=PFILES;
  int pr_mzintegral=PMZINT;
  double mz_int_arr[NPOINTS+1][2]={0};
  int mz_int_counter=0;
  countername=(char *)malloc(sizeof(char)*2);
  float radius=A;
  float  veloc=V;


  f0 = ( double * ) malloc ( m * sizeof ( double ) );
  f1 = ( double * ) malloc ( m * sizeof ( double ) );
  f2 = ( double * ) malloc ( m * sizeof ( double ) );
  f3 = ( double * ) malloc ( m * sizeof ( double ) );
  u0 = ( double * ) malloc ( m * sizeof ( double ) );
  u1 = ( double * ) malloc ( m * sizeof ( double ) );
  u2 = ( double * ) malloc ( m * sizeof ( double ) );
  u3 = ( double * ) malloc ( m * sizeof ( double ) );



  dt = ( tspan[1] - tspan[0] ) / ( double ) ( n );
  points=n/NPOINTS;
  tpoints=points;
  j = 0;
  t[0] = tspan[0];

  if(prfiles==1){
    sprintf(countername,"%hu",(unsigned short)0);
    countername[1]='\0';
    // strcat(filename,"./data/vort");
    strcpy(filename,"/home/polychronou/Documents/data/vort\0");
    strcat(filename, countername);
    strcat(filename, ".txt");
    f=fopen(filename,"w+");

    if(f==NULL){
      errnum = errno;
      fprintf(stderr, "Value of errno: %d\n", errno);
      perror("Error printed by perror open");
      fprintf(stderr, "Error opening file: %s\n", strerror( errnum ));
      exit(1);
    }
  }

  for ( i = 0; i < m; i++ )
  {
    y[i+j*m] = y0[i];
    y[i+m]=y0[i];
    if (prfiles) {fprintf(f, "%lf ", y0[i]);}
    f0[i]=0;f1[i]=0;f2[i]=0;f3[i]=0;
    u0[i]=0;u1[i]=0;u2[i]=0;u3[i]=0;
  }
  time_point=tpoints*1e-3+tspan[0];
  if (prfiles==1){
    fprintf(f,"\n");
    fprintf(f,"Ex=%lf  E_an=%lf  E_dm=%lf  sum=%lf",EXenergy((y0),r,z,dr,dz),
    ANenergy((y+m),r,z,dr,dz),DMenergy(y+m,r,z,dr,dz),EXenergy((y0),r,z,dr,dz)+
    ANenergy((y+m),r,z,dr,dz)+DMenergy(y+m,r,z,dr,dz));
    fprintf(f,"\n");
    fprintf(f,"t=%1.2lf\n",0.00);
    filecounter++;
    fclose(f);
  }
  if (pr_mzintegral==1){
    mz_int_arr[mz_int_counter][0]=mz_integral(y+m,r,z,dr,dz);
    mz_int_arr[mz_int_counter][1]=0.00;
    mz_int_counter++;
  }

  for ( j = 0; j < n; j++ )
  {
    t0 = t[j];
    norm=0;
    for ( i = 0; i < m; i++ )
    {
      // norm+=pow((y[i+j*m*0]-y[i+(j*0+1)*m]),2);
      y[i+j*m*0]=y[i+(j*0+1)*m];
      u0[i] = y[i+j*m*0];

    }
    if(j==0){
      ;
    }else{
      preven=en;
    }
    // norm=sqrt(norm);
    // if(norm<1e-1 && j>0){
    //   printf("norm: %lf\n",norm);
    //   break;
    // }
    // if (j>28){
    //   ;
    // }
    dydt ( t0, u0, f0, r,z, dr,dz );

    t1 = t0 + dt / 2.0;
    for ( i = 0; i < m; i++ )
    {
      u1[i] = u0[i] + dt * f0[i] / 2.0;
    }
    dydt ( t1, u1, f1, r ,z, dr,dz );

    t2 = t0 + dt / 2.0;
    for ( i = 0; i < m; i++ )
    {
      u2[i] = u0[i] + dt * f1[i] / 2.0;
    }
    dydt ( t2, u2, f2, r,z ,dr,dz );

    t3 = t0 + dt;
    for ( i = 0; i < m; i++ )
    {
      u3[i] = u0[i] + dt * f2[i];
    }
    dydt ( t3, u3, f3 ,r,z,dr,dz);

    t[j+1] = t[j] + dt;
    for ( i = 0; i < m; i++ )
    {
       y[i+(j*0+1)*m] = u0[i] + dt * ( f0[i] + 2.0 * f1[i] + 2.0 * f2[i] + f3[i] ) / 6.0;
    }
    norm=0;
    if(j==tpoints ){
      if(prfiles==1){
        if(filecounter>=10){
          countername=(char*)realloc(countername,sizeof(char)*3);
          sprintf(countername,"%hu",(unsigned short)filecounter);
          countername[2]='\0';
        }else{
        sprintf(countername,"%hu",(unsigned short)filecounter);
        countername[1]='\0';
        }
        // strcat(filename,"./data/vort");
        strcpy(filename,"/home/polychronou/Documents/data/vort\0");
        strcat(filename, countername);
        strcat(filename, ".txt\0");
        f=fopen(filename,"w+");

        if(f==NULL){
          errnum = errno;
          fprintf(stderr, "Value of errno: %d\n", errno);
          perror("Error printed by perror open");
          fprintf(stderr, "Error opening file: %s\n", strerror( errnum ));
          exit(1);
        }
      }


    }

    for(i=0; i<m; i=i+3){
      norm=sqrt(pow(y[i+(j*0+1)*m],2)+pow(y[i+1+(j*0+1)*m],2)+pow(y[i+2+(j*0+1)*m],2));
      y[i+(j*0+1)*m]/=norm;
      y[i+1+(j*0+1)*m]/=norm;
      y[i+2+(j*0+1)*m]/=norm;
      if(j==tpoints && prfiles==1){
        fprintf(f, "%lf ", y[i+0+(j*0+1)*m]);
        fprintf(f, "%lf ", y[i+1+(j*0+1)*m]);
        fprintf(f, "%lf ", y[i+2+(j*0+1)*m]);
      }
    }
    e_ex=EXenergy(y+m,r,z,dr,dz);
    e_dm=DMenergy(y+m,r,z,dr,dz);
    e_an=ANenergy(y+m,r,z,dr,dz);
    en=e_ex+e_dm+e_an;
    if(j==tpoints ){
      time_point=tpoints*1e-3+tspan[0];
      tpoints+=points;

      if (prfiles==1){
        fprintf(f,"\n");
        fprintf(f,"Ex=%lf  E_an=%lf  E_dm=%lf  sum=%lf",e_ex,e_an,e_dm,e_ex+e_an+e_dm);
        fprintf(f,"\n");
        fprintf(f,"t=%1.2lf\n",time_point);
        filecounter++;
        fclose(f);
      }


      if(fabs((en-preven)/preven)<1e-12){
        printf("No energy variation en=%lf , preven=%lf\n",en,preven);
        break;
      }

    if (pr_mzintegral==1){
      mz_int_arr[mz_int_counter][0]=mz_integral(y+m,r,z,dr,dz);
      mz_int_arr[mz_int_counter][1]=time_point;
      mz_int_counter++;
    }
  }
  }

  if (pr_mzintegral==1){
      f =fopen("/home/polychronou/Documents/data/mz_integral.txt","w+");
      if(f==NULL){
        errnum = errno;
        fprintf(stderr, "Value of errno: %d\n", errno);
        perror("Error printed by perror open");
        fprintf(stderr, "Error opening file: %s\n", strerror( errnum ));
        exit(1);
      }
      for (i=0;i<NPOINTS+1;i++){
        if(i>0){
          fprintf(f, " ");
        }
        fprintf(f,"%lf %lf",mz_int_arr[i][0],mz_int_arr[i][1]);
      }
        fprintf(f,"\n");
        fprintf(f,"radius:%1.3f velocity:%1.3f ",radius,veloc);
      fclose(f);
  }

/*
  Free memory.
*/
  free ( f0 );
  free ( f1 );
  free ( f2 );
  free ( f3 );
  free ( u0 );
  free ( u1 );
  free ( u2 );
  free ( u3 );
  free (countername);
  return;
}
/******************************************************************************/

void timestamp ( )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 September 2003

  Author:

    John Burkardt
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
