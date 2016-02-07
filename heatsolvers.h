typedef struct _simulation
{
 double ***grid;

 double d_x;
 double d_y;
 double d_z;
 
 double l_x;
 double l_y;
 double l_z;

 int Nx;
 int Ny;
 int Nz;

 double alpha;
 
 double (*S)(double x, double y, double z);

 double low_x;
 double high_x;
 double low_y;
 double high_y;
 double low_z;
 double high_z;

 short periodic; //boolean: 1 means it is periodic, 0 means Dirichlet
} simulation;

void gauss_init_simulation(simulation *A, double Amp, double x0, double y0, double z0, double sigx, double sigy, double sigz);

void print_sim(simulation *A, int x, FILE *fp);

void cranknicholson(simulation *A, int nt, double dt);
void ftcs(simulation *A, int nt, double dt);
void adi(simulation *A, int nt, double dt);
void rref(double **A, int m, int n);
