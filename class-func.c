#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<math.h>

/* a function to get a random double between min and max*/
double get_double(double min, double max, unsigned int *seed)
{
 double val = min;
 double scope = (max-min);
 double inc = (scope/RAND_MAX);
 val += (rand_r(seed)*inc);
 return val;
}

double get_time(struct timespec start, struct timespec stop)
{
 return(
	((double)stop.tv_sec + (double)stop.tv_nsec/1000000000) -
	((double)start.tv_sec + (double)start.tv_nsec/1000000000)
	);
}

void print_data(double *data, int size, FILE *fp)
{
 int i;

 for(i = 0; i < size; ++i)
  fprintf(fp, "%.20lf\t", data[i]);
 fprintf(fp, "\n");

 return;
}


