#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include"class-func.h"
#include"heatsolvers.h"
#include"nrutil.h"
#define NUMTESTS 5
#define START 3

int main()
{
 struct timespec start, stop; 
 int n, i;
 double data[NUMTESTS][4];

 FILE *fp = fopen("results", "w");
 simulation *A = malloc(sizeof(simulation));
 simulation *B = malloc(sizeof(simulation));
 simulation *C = malloc(sizeof(simulation));

 for(i = 0; i < NUMTESTS; ++i)
 {
  n = START+i;
  A->l_x = A->l_y = A->l_z = n;
  A->d_x = A->d_y = A->d_z = 1;
  A->alpha = 0.1;
  A->periodic = 0;
  A->high_x = A->low_x = A->high_z = A->low_z = A->high_y = A->low_y = 0;
  A->S = NULL;
  gauss_init_simulation(A,1,n/2,n/2,n/2,n/3,n/3,n/3);

  B->l_x = B->l_y = B->l_z = n;
  B->d_x = B->d_y = B->d_z = 1;
  B->alpha = 0.1;
  B->periodic = 0;
  B->high_x = B->low_x = B->high_z = B->low_z = B->high_y = B->low_y = 0;
  B->S = NULL;
  gauss_init_simulation(B,1,n/2,n/2,n/2,n/3,n/3,n/3);

  C->l_x = C->l_y = C->l_z = n;
  C->d_x = C->d_y = C->d_z = 1;
  C->alpha = 0.1;
  C->periodic = 0;
  C->high_x = C->low_x = C->high_z = C->low_z = C->high_y = C->low_y = 0;
  C->S = NULL;
  gauss_init_simulation(C,1,n/2,n/2,n/2,n/3,n/3,n/3);

  clock_gettime(CLOCK_MONOTONIC, &start);
  ftcs(A, 100, 0.1);
  clock_gettime(CLOCK_MONOTONIC, &stop);
  data[i][1] = get_time(start, stop);

  clock_gettime(CLOCK_MONOTONIC, &start);
  cranknicholson(B, 100, 0.1);
  clock_gettime(CLOCK_MONOTONIC, &stop);
  data[i][2] = get_time(start, stop);

  clock_gettime(CLOCK_MONOTONIC, &start);
  adi(C, 100, 0.1);
  clock_gettime(CLOCK_MONOTONIC, &stop);
  data[i][3] = get_time(start, stop);

  data[i][0] = n;
 
  free(A->grid);
  free(B->grid);
  free(C->grid);
 }

 for(i = 0; i < NUMTESTS; ++i)
  print_data(data[i], 4, fp);

 fclose(fp);

 free(A);
 free(B);
 free(C);
 return 0;
}
