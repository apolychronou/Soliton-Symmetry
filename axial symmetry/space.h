float deriv2r(double *m,int i,double dr);
float deriv1r(double *m,int i, double dr);
float deriv2z(double *m,int i,double dz);
float deriv1r(double *m,int i,double dz);
void laplace(double *m,double *r,double *z,double *laplacian,int ri,int zi,double dr,double dz);
void skyrmion_init_values(double *m0,double *r, double *z);
void LLequation(double t, double *m,double *dm,double *r,double *z,double dr,double dz);
double EXenergy(double *m, double *r,double *z, double dr,double dz);
double ANenergy(double *m, double *r,double *z, double dr,double dz);
double DMenergy(double *m, double *r,double *z, double dr,double dz);