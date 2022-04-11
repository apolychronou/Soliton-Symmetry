#define N 300 //r grid
#define NZ 200 //z grid
#define DL 0.6 //d-m interaction coef
#define DP 0.01//dumping coef
#define AN 1  //anisotropy coef

int magnIndex(int i, int j);
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
void DMenergy_Density(double *m, double *r,double *z, double dr,double dz,double a,double phi);
void center_angle(double *y,double *r,double *z, double dr, double dz);
void DMenergy_zDensity(double *m, double *r,double *z, double dr,double dz,double a,double phi);
void skyrmion_init_values(double *m0,double *r,double *z);
void boundary_conditions(double *m);
