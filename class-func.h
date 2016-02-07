#include<time.h>
#include<stdio.h>

double get_double(double min, double max, unsigned int *seed);
double get_time(struct timespec start, struct timespec stop);
void print_data(double *data, int size, FILE *fp);
