# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <time.h>
#include <string.h>
#include <errno.h>
# include "space.h"
# include "rk4.h"

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
  double norm=0;
  int points,tpoints;
  int filecounter=0;
  char *countername;
  char filename[17]="./data/vort";
  FILE *f;
  extern int errno;
  errno=0;
  int errnum;

  countername=(char *)malloc(sizeof(char)*1);


  f0 = ( double * ) malloc ( m * sizeof ( double ) );
  f1 = ( double * ) malloc ( m * sizeof ( double ) );
  f2 = ( double * ) malloc ( m * sizeof ( double ) );
  f3 = ( double * ) malloc ( m * sizeof ( double ) );
  u0 = ( double * ) malloc ( m * sizeof ( double ) );
  u1 = ( double * ) malloc ( m * sizeof ( double ) );
  u2 = ( double * ) malloc ( m * sizeof ( double ) );
  u3 = ( double * ) malloc ( m * sizeof ( double ) );

  dt = ( tspan[1] - tspan[0] ) / ( double ) ( n );
  points=n/10;
  tpoints=points;
  j = 0;
  t[0] = tspan[0];

  sprintf(countername,"%d",filecounter);
  // strcat(filename,"./data/vort");
  strcpy(filename,"./data/vort\0");
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


  for ( i = 0; i < m; i++ )
  {
    y[i+j*m] = y0[i];
    y[i+m]=y0[i];
    fprintf(f, "%lf ", y0[i]);
  }
  fprintf(f,"\n");
  fprintf(f,"Ex=%lf  E_an=%lf  E_dm=%lf  sum=%lf",EXenergy((y0),r,z,dr,dz),
  ANenergy((y+m),r,z,dr,dz),DMenergy(y+m,r,z,dr,dz),EXenergy((y0),r,z,dr,dz)+
  3*ANenergy((y+m),r,z,dr,dz)+2*DMenergy(y+m,r,z,dr,dz));
  filecounter++;
  fclose(f);

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
    if(j==tpoints){
      if(filecounter>10){
        countername=(char*)realloc(countername,sizeof(char)*2);
      }
      sprintf(countername,"%d",filecounter);
      // strcat(filename,"./data/vort");
      strcpy(filename,"./data/vort\0");
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

    for(i=0; i<m; i=i+3){
      norm=sqrt(pow(y[i+(j*0+1)*m],2)+pow(y[i+1+(j*0+1)*m],2)+pow(y[i+2+(j*0+1)*m],2));
      y[i+(j*0+1)*m]/=norm;
      y[i+1+(j*0+1)*m]/=norm;
      y[i+2+(j*0+1)*m]/=norm;
      if(j==tpoints){
        fprintf(f, "%lf ", y[i+0+(j*0+1)*m]);
        fprintf(f, "%lf ", y[i+1+(j*0+1)*m]);
        fprintf(f, "%lf ", y[i+2+(j*0+1)*m]);
      }
    }

    if(j==tpoints){
      fprintf(f,"\n");
      fprintf(f,"Ex=%lf  E_an=%lf  E_dm=%lf  sum=%lf",EXenergy((y+m),r,z,dr,dz),
      ANenergy((y+m),r,z,dr,dz),DMenergy(y+m,r,z,dr,dz),EXenergy((y+m),r,z,dr,dz)+
      3*ANenergy((y+m),r,z,dr,dz)+2*DMenergy(y+m,r,z,dr,dz));
      filecounter++;
      fclose(f);
      tpoints+=points;
    }


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
