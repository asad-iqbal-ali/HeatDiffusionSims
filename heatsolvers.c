#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<string.h>
#include"heatsolvers.h"
#include"nrutil.h"

#define SWAP(a,b) {double temp=(a);(a)=(b);(b)=temp;}
#define DEBUG 0

void gauss_init_simulation(simulation *A, double Amp, double x0, double y0, double z0, double sigx, double sigy, double sigz)
{
 int i, j, k;
 A->Nx = A->l_x/A->d_x;
 A->Ny = A->l_y/A->d_y;
 A->Nz = A->l_z/A->d_z;
 
 double sigsqx = 2*sigx*sigx;
 double sigsqy = 2*sigy*sigy;
 double sigsqz = 2*sigz*sigz;

 double x, y, z; 


 /*malloc out the actual grid*/
 A->grid = (double ***)malloc(sizeof(double **)*A->Nx);
 assert(A->grid != NULL);
 double **data = (double **)malloc(sizeof(double *)*A->Nx*A->Ny);
 assert(data != NULL);
 double *bigdata = malloc(sizeof(double)*A->Nx*A->Ny*A->Nz); 
 assert(bigdata != NULL);
 
 for(i = 0; i < A->Nx; ++i)
 {
  A->grid[i] = data+(i*A->Ny);
  for(j = 0; j < A->Ny; ++j)
   (data+(i*A->Ny))[j] = bigdata+(A->Nz*(j+(i*A->Ny)));
 }

 /*insert values into the grid by gauss+noise*/

 for(i = 0; i < A->Nx; ++i)
 {
  x = i*A->d_x;
  for(j = 0; j < A->Ny; ++j)
  {
   y = j*A->d_y;
   for(k = 0; k < A->Nz; ++k)
   {
    z = k*A->d_z;
    A->grid[i][j][k] = Amp*exp(-
	(
	 (x-x0)*(x-x0)/sigsqx +
	 (y-y0)*(y-y0)/sigsqy +
	 (z-z0)*(z-z0)/sigsqz
	))
	+(A->S == NULL? 0.0:A->S(x,y,z));
   }
  }
 } 
}

/*a function to print a 2D slice of the simulation at (confusingly) z-value x*/
void print_sim(simulation *A, int x, FILE *fp)
{
 int j, k;
 assert(x > 0);
 assert(x < (A->l_x/A->d_x));
 assert(fp != NULL);

 for(j = 0; j < A->Ny; ++j)
 {
  for(k = 0; k < A->Nz; ++k)
   fprintf(fp, "%lf ", A->grid[x][j][k]);
  fprintf(fp, "\n");
 }
}



void cranknicholson(simulation *sim, int nt, double dt)
{
 int i,j,k;
 int x,y,z;
 int t;
 int size = sim->Nx*sim->Nz*sim->Ny;

 double **A,**augA;
 double *Adata, *augAdata, **bdata, *bbigdata;
 double up, down, front, back, left, right;
 double C = (sim->alpha*dt)/(2*(sim->d_x*sim->d_x));

 
 /*initialize A as a sizexsize matrix of zeroes*/
 A = (double **)malloc(sizeof(double *)*size);
 assert(A != NULL);
 Adata = malloc(sizeof(double)*size*size); 
 assert(Adata != NULL);
 
 for(i = 0; i < size; ++i)
 {
  A[i] = Adata+(i*size);
  memset(A[i], 0.0, sizeof(double)*size);
 }


 /*augA is the augmented matrix of A and b that we will rref*/ 
 augA = (double **)malloc(sizeof(double *)*size);
 assert(augA != NULL);
 augAdata = malloc(sizeof(double)*size*(size+1)); 
 assert(augAdata != NULL);
 
 for(i = 0; i < size; ++i)
 {
  augA[i] = augAdata+(i*(size+1));
 }

 /*set non-zero values in A*/
 for(i = 0; i < size; ++i)
 {
  /*The analagous x,y,and z values for a given spot in a row*/
  x = i/(sim->Ny*sim->Nz);
  y = (i/sim->Nz)%sim->Ny;
  z = i%sim->Nz;

  /*diagonal values; coefficient for T(i,j,k) for the next timestep*/
  A[i][i] = (1+(6.0*C));

  /*place "-C" values in their appropriate spots,
   accounting for border values. If the value is on a border,
   include wraparound if the boundaries are periodic. If they are not,
   then border values are a known quantity and therefore not
   included in A*/
 
  
  if(!x)
  {
   if(sim->periodic)
    A[i][i+((sim->Nx-1)*sim->Ny*sim->Nz)] = -C;
  }
  else A[i][i - (sim->Ny*sim->Nz)] = -C;
  if(x == sim->Nx-1)
  {
   if(sim->periodic)
    A[i][i - ((sim->Nx-1)*sim->Ny*sim->Nz)] = -C;
  }
  else A[i][i + (sim->Ny*sim->Nz)] = -C;
 
  if(!y)
  {
   if(sim->periodic)
    A[i][i+((sim->Ny-1)*sim->Nz)] = -C;
  }
  else A[i][i - sim->Nz] = -C;
  if(y == sim->Ny-1)
  {
   if(sim->periodic)
    A[i][i - ((sim->Ny-1)*sim->Nz)] = -C;
  }
  else A[i][i + sim->Nz] = -C;


  if(!z)
  {
   if(sim->periodic)
    A[i][i+sim->Nz-1] = -C;
  }
  else A[i][i - 1] = -C;
  if(z == sim->Nz-1)
  {
   if(sim->periodic)
    A[i][i - (sim->Nz-1)] = -C;
  }
  else A[i][i+1] = -C;
 
 
 }

/*The actual simulation*/

 for(t = 0; t < nt; ++t)
 {
  /*First, fill in the values in augA from A*/
  for(i = 0; i < size; ++i)
   memcpy(augA[i],A[i],sizeof(double)*size);
  /*Next, fill in the values for 'b'*/
  for(i = 0; i < sim->Nx; ++i)
  {
   for(j = 0; j < sim->Ny; ++j)
   {
    for(k = 0; k < sim->Nz; ++k)
    {
     /*figure out the edge values.*/
     /*If we're at i=0, first check if the grid is periodic*/
     left = (i? sim->grid[i-1][j][k] : 
     /*if so, use the value for i = Nx-1. Otherwise use the specified boundary value.
	in this situation, the boundary term would also be included from the next time step,
	since it would be known, so we include twice the known boundary value*/
		(sim->periodic? sim->grid[sim->Nx - 1][j][k] : 2.0*sim->low_x));
     /*Repeat similarly structured checks for j and k*/
     right = (i< (sim->Nx - 1)? sim->grid[i+1][j][k] : 
		(sim->periodic? sim->grid[0][j][k] : 2.0*sim->high_x));
     up = (j? sim->grid[i][j-1][k] : 
		(sim->periodic? sim->grid[i][sim->Ny - 1][k] : 2.0*sim->low_y));
     down = (j< (sim->Ny - 1)? sim->grid[i][j+1][k] : 
		(sim->periodic? sim->grid[i][0][k] : 2.0*sim->high_y));
     front = (k? sim->grid[i][j][k-1] :
		(sim->periodic? sim->grid[i][j][sim->Nz - 1] :2.0*sim->low_z));
     back = (k< (sim->Nz - 1)? sim->grid[i][j][k+1] : 
		(sim->periodic? sim->grid[i][j][0] : 2.0*sim->high_z));
     /*Fill in the value in 'b' for AugA as per the Crank-Nicolson equation*/
     augA[(i*sim->Ny*sim->Nz)+(j*(sim->Nz))+k][size] = 
	((1.0-(6.0*C))*sim->grid[i][j][k]) +
	(C*(up+down+left+right+front+back)) +
	((sim->S == NULL)? 0.0 : sim->S(i*sim->d_x,j*sim->d_y,k*sim->d_z));
    }
   }	
  }
  rref(augA,size,size+1);
  
  /*copy over the new T-values into the grid*/
  for(i = 0; i < size; ++i)
   sim->grid[0][0][i] = augA[i][size];

 }

 free(Adata);
 free(A);
 free(augAdata);
 free(augA);

}

void ftcs(simulation *A, int nt, double dt)
{
 int t, i, j, k;
 double C = dt*A->alpha;
 double ***tmp;
 double dxsq = A->d_x*A->d_x;
 double dysq = A->d_y*A->d_y;
 double dzsq = A->d_z*A->d_z;
 double left, right, up, down, front, back;

 /*Stability conditions: assume that this has to hold true
   for all dx, dy, and dz*/
 assert((C/dxsq) < 0.125);
 assert((C/dysq) < 0.125);
 assert((C/dzsq) < 0.125);

 /*create the second grid, which will be used to calculate through the timesteps*/
 double ***newgrid = (double ***)malloc(sizeof(double **)*A->Nx);
 assert(newgrid != NULL);
 double **data = (double **)malloc(sizeof(double *)*A->Nx*A->Ny);
 assert(data != NULL);
 double *bigdata = malloc(sizeof(double)*A->Nx*A->Ny*A->Nz); 
 assert(bigdata != NULL);
 
 for(i = 0; i < A->Nx; ++i)
 {
  newgrid[i] = data+(i*A->Ny);
  for(j = 0; j < A->Ny; ++j)
   (data+(i*A->Ny))[j] = bigdata+(A->Nz*(j+(i*A->Ny)));
 }

 for(t = 0; t < nt; ++t)
 {
  for(i = 0; i < A->Nx; ++i)
   for(j = 0; j < A->Ny; ++j)
    for(k = 0; k < A->Nz; ++k)
    {
     /*figure out the edge values for the FTCS equation.*/

     /*If we're at i=0, first check if the grid is periodic*/
     left = (i? A->grid[i-1][j][k] : 
     /*if so, use the value for i = Nx-1. Otherwise use the specified value*/
		(A->periodic? A->grid[A->Nx - 1][j][k] : A->low_x));
     /*Repeat similarly structured checks for j and k*/
     right = (i< (A->Nx - 1)? A->grid[i+1][j][k] : 
		(A->periodic? A->grid[0][j][k] : A->high_x));
     up = (j? A->grid[i][j-1][k] : 
		(A->periodic? A->grid[i][A->Ny - 1][k] : A->low_y));
     down = (j< (A->Ny - 1)? A->grid[i][j+1][k] : 
		(A->periodic? A->grid[i][0][k] : A->high_y));
     front = (k? A->grid[i][j][k-1] :
		(A->periodic? A->grid[i][j][A->Nz - 1] : A->low_z));
     back = (k< (A->Nz - 1)? A->grid[i][j][k+1] : 
		(A->periodic? A->grid[i][j][0] : A->high_z));
     /*Do the actual calculation*/
     newgrid[i][j][k] = A->grid[i][j][k] + C*
	(
	 (left - 2*A->grid[i][j][k] + right)/dxsq
	+(up - 2*A->grid[i][j][k] + down)/dysq
	+(front - 2*A->grid[i][j][k] + back)/dzsq
	)
	+(A->S == NULL? 0.0:A->S(i*A->d_x, j*A->d_y, k*A->d_z));
    }
  tmp = A->grid;
  A->grid = newgrid;
  newgrid = tmp;
 }

 free(newgrid[0][0]);
 free(newgrid[0]);
 free(newgrid);

}


void adi(simulation *sim, int nt, double dt)
{
 int i,j,k;
 int x,y,z;
 int t;
 int size = sim->Nx*sim->Nz*sim->Ny;

 double **A[3],**augA;
 double *Adata[3], *augAdata, **bdata, *bbigdata;
 double up, down, front, back, left, right;
 double C = (sim->alpha*dt)/(3.0*(sim->d_x*sim->d_x));

 
 /*initialize the three A's as a sizexsize matrix of zeroes
  One A for each partial timestep*/
 for(t = 0; t < 3; ++t)
 {
  A[t] = (double **)malloc(sizeof(double *)*size);
  assert(A[t] != NULL);
  Adata[t] = malloc(sizeof(double)*size*size); 
  assert(Adata[t] != NULL);
  
  for(i = 0; i < size; ++i)
  {
   A[t][i] = Adata[t]+(i*size);
   memset(A[t][i], 0.0, sizeof(double)*size);
  }
 }

 /*augA is the augmented matrix of A and b that we will rref*/ 
 augA = (double **)malloc(sizeof(double *)*size);
 assert(augA != NULL);
 augAdata = malloc(sizeof(double)*size*(size+1)); 
 assert(augAdata != NULL);
 
 for(i = 0; i < size; ++i)
 {
  augA[i] = augAdata+(i*(size+1));
 }

 /*set non-zero values in A's*/
 for(i = 0; i < size; ++i)
 {
  /*The analagous x,y,and z values for a given spot in a row*/
  x = i/(sim->Ny*sim->Nz);
  y = (i/sim->Nz)%sim->Ny;
  z = i%sim->Nz;

  /*diagonal values; coefficient for T(i,j,k) for the next timestep*/
  A[0][i][i] = A[1][i][i] = A[2][i][i] = 1.0+(2.0*C);

  /*place "-C" values in their appropriate spots,
   accounting for border values. If the value is on a border,
   include wraparound if the boundaries are periodic. If they are not,
   then border values are a known quantity and therefore not
   included in A. Each A only deals with the neighbors in the dimension
   for which it will calculate values. These are x, y, and z, for A0, A1,
   and A2, respectively*/
 
  if(!x)
  {
   if(sim->periodic)
    A[0][i][i+((sim->Nx-1)*sim->Ny*sim->Nz)] = -C;
  }
  else A[0][i][i - (sim->Ny*sim->Nz)] = -C;
  if(x == sim->Nx-1)
  {
   if(sim->periodic)
    A[0][i][i - ((sim->Nx-1)*sim->Ny*sim->Nz)] = -C;
  }
  else A[0][i][i + (sim->Ny*sim->Nz)] = -C;
 
  if(!y)
  {
   if(sim->periodic)
    A[1][i][i+((sim->Ny-1)*sim->Nz)] = -C;
  }
  else A[1][i][i - sim->Nz] = -C;
  if(y == sim->Ny-1)
  {
   if(sim->periodic)
    A[1][i][i - ((sim->Ny-1)*sim->Nz)] = -C;
  }
  else A[1][i][i + sim->Nz] = -C;


  if(!z)
  {
   if(sim->periodic)
    A[2][i][i+sim->Nz-1] = -C;
  }
  else A[2][i][i - 1] = -C;
  if(z == sim->Nz-1)
  {
   if(sim->periodic)
    A[2][i][i - (sim->Nz-1)] = -C;
  }
  else A[2][i][i+1] = -C;
 
 
 }
/*The actual simulation*/

 for(t = 0; t < nt; ++t)
 {
  /*First, fill in the values in augA from A*/
  for(i = 0; i < size; ++i)
   memcpy(augA[i],A[0][i],sizeof(double)*size);
  /*Next, fill in the values for 'b'*/
  for(i = 0; i < sim->Nx; ++i)
  {
   for(j = 0; j < sim->Ny; ++j)
   {
    for(k = 0; k < sim->Nz; ++k)
    {
     /*This works the same as the other two functions, except that in the dimension
	being calculated (in this instance, x) we only need to include neighbors for that
	dimension in the event they are known, which is only if they are border values
	and the simulation is not periodic. Outside of this case, use 0.0.*/
     left = (i||sim->periodic? 0.0 : sim->low_x);
     right = (i< (sim->Nx - 1)||sim->periodic? 0.0 : sim->high_x);
     up = (j? sim->grid[i][j-1][k] : 
		(sim->periodic? sim->grid[i][sim->Ny - 1][k] : 2.0*sim->low_y));
     down = (j< (sim->Ny - 1)? sim->grid[i][j+1][k] : 
		(sim->periodic? sim->grid[i][0][k] : 2.0*sim->high_y));
     front = (k? sim->grid[i][j][k-1] :
		(sim->periodic? sim->grid[i][j][sim->Nz - 1] :2.0*sim->low_z));
     back = (k< (sim->Nz - 1)? sim->grid[i][j][k+1] : 
		(sim->periodic? sim->grid[i][j][0] : 2.0*sim->high_z));
     augA[(i*sim->Ny*sim->Nz)+(j*(sim->Nz))+k][size] = 
	((1.0-(4.0*C))*sim->grid[i][j][k]) +
	(C*(up+down+left+right+front+back)) +
	((sim->S == NULL)? 0.0 : sim->S(i*sim->d_x,j*sim->d_y,k*sim->d_z)/3.0);
    }
   }	
  }
  rref(augA,size,size+1);
  for(i = 0; i < size; ++i)
   sim->grid[0][0][i] = augA[i][size];

  /*repeat the steps for the y dimension*/

  for(i = 0; i < size; ++i)
   memcpy(augA[i],A[1][i],sizeof(double)*size);
  for(i = 0; i < sim->Nx; ++i)
  {
   for(j = 0; j < sim->Ny; ++j)
   {
    for(k = 0; k < sim->Nz; ++k)
    {
     left = (i? sim->grid[i-1][j][k] : 
		(sim->periodic? sim->grid[sim->Nx - 1][j][k] : 2.0*sim->low_x));
     right = (i< (sim->Nx - 1)? sim->grid[i+1][j][k] : 
		(sim->periodic? sim->grid[0][j][k] : 2.0*sim->high_x));
     up = (j||sim->periodic? 0.0: sim->low_y);
     down = (j< (sim->Ny - 1)||sim->periodic? 0.0 : sim->high_y);
     front = (k? sim->grid[i][j][k-1] :
		(sim->periodic? sim->grid[i][j][sim->Nz - 1] :2.0*sim->low_z));
     back = (k< (sim->Nz - 1)? sim->grid[i][j][k+1] : 
		(sim->periodic? sim->grid[i][j][0] : 2.0*sim->high_z));
     augA[(i*sim->Ny*sim->Nz)+(j*(sim->Nz))+k][size] = 
	((1.0-(4.0*C))*sim->grid[i][j][k]) +
	(C*(up+down+left+right+front+back)) +
	((sim->S == NULL)? 0.0 : sim->S(i*sim->d_x,j*sim->d_y,k*sim->d_z)/3.0);
    }
   }	
  }
  rref(augA,size,size+1);
  for(i = 0; i < size; ++i)
   sim->grid[0][0][i] = augA[i][size];
  
  /*repeat again for the z dimension*/
  for(i = 0; i < size; ++i)
   memcpy(augA[i],A[2][i],sizeof(double)*size);
  for(i = 0; i < sim->Nx; ++i)
  {
   for(j = 0; j < sim->Ny; ++j)
   {
    for(k = 0; k < sim->Nz; ++k)
    {
     left = (i? sim->grid[i-1][j][k] : 
		(sim->periodic? sim->grid[sim->Nx - 1][j][k] : 2.0*sim->low_x));
     right = (i< (sim->Nx - 1)? sim->grid[i+1][j][k] : 
		(sim->periodic? sim->grid[0][j][k] : 2.0*sim->high_x));
     up = (j? sim->grid[i][j-1][k] : 
		(sim->periodic? sim->grid[i][sim->Ny - 1][k] : 2.0*sim->low_y));
     down = (j< (sim->Ny - 1)? sim->grid[i][j+1][k] : 
		(sim->periodic? sim->grid[i][0][k] : 2.0*sim->high_y));
     front = (k ||sim->periodic? 0.0 :sim->low_z);
     back = (k< (sim->Nz - 1)||sim->periodic? 0.0: sim->high_z);
     augA[(i*sim->Ny*sim->Nz)+(j*(sim->Nz))+k][size] = 
	((1.0-(4.0*C))*sim->grid[i][j][k]) +
	(C*(up+down+left+right+front+back)) +
	((sim->S == NULL)? 0.0 : sim->S(i*sim->d_x,j*sim->d_y,k*sim->d_z)/3.0);
    }
   }	
  }
  rref(augA,size,size+1);
  for(i = 0; i < size; ++i)
   sim->grid[0][0][i] = augA[i][size];
 }

 for(i = 0; i < 3; ++i)
 {
  free(Adata[i]);
  free(A[i]);
 }
  free(augAdata);
  free(augA);
}

void rref(double **A, int m, int n)
{
  double* tmp=dvector(0,n-1);
  double tol=0.000000000000001;
  double p;
  long i,j,k,l;
  i=0; 
  j=0;
  while ((i < m) && (j< n)){
    /* Initialize pivot */
    p=0.; k=i;
    /* Extract best pivot, and it's row index 
       (even non-zero entries.  Called partial pivoting -- improves numerical stability)
    */
    for(l=i;l<m;l++){
      if ( p < fabs( A[l][j] ) ){
	p=A[l][j];
	k=l;
      }
    }
    if (fabs(p)<=tol){
      for(l=i;l<m;l++){
	A[l][j]=0.;
      }
      j++;
      continue;
    }
    if(k!=i) /* Swap ith and kth rows if needed */
      for(l=j;l<n;l++) SWAP( A[i][l] , A[k][l] );
    
    /* Divide pivot row by pivot element */
    p=A[i][j];
    for (l=j;l<n;l++) A[i][l] = A[i][l]/p;
    
    /* Subtract multiples of the pivot row from all
       other rows  */
    for (k=0;k<m;k++){
      if (k!=i){
	for (l=j;l<n;l++)
	  tmp[l]=A[k][j]*A[i][l];
	for (l=j;l<n;l++)
	  A[k][l]=A[k][l]-tmp[l];
      }
    }
    if (DEBUG)
      mprint(A,m,n,"A intermediate");
    i++; 
    j++;
  }
}

