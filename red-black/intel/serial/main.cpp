#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <omp.h>
#include <time.h>
#include "core.h"
#include "util.h"

#include <unistd.h>

#define L 500

#define PI 3.14159265359

/*
static double vertical_boundary(double x)
{
	UNUSED(x);
	return 0;
}

static double horizontal_boundary(double x)
{
	UNUSED(x);
	return 0;
}
*/

double timeInSeconds(struct timespec *t);

int main()
{
	struct timespec time1, time2;
	// omp_set_num_threads(8);
	clock_gettime(CLOCK_MONOTONIC, &time1);

	double *T = make_grid(L);
	double *b = make_grid(L);

	GRID_GET(b, L/4, L/4, L) = 1;
	//GRID_GET(b, 3*L/4, L/4, L) = 1;
	//GRID_GET(b, L/4, 3*L/4, L) = 1;
	GRID_GET(b, 3*L/4, 3*L/4, L) = 1;

	set_boundary(T, L, std::cos, std::sin, std::cos, std::sin);
	iterate_red_black(T, b, L, 10000, 0.5);
	dump_grid(T, L, "out.csv");

	clock_gettime(CLOCK_MONOTONIC, &time2);
    double testTime = timeInSeconds(&time2)-timeInSeconds(&time1);
    printf("Test time: %f\n", testTime);
}

double timeInSeconds(struct timespec *t){
    // a timespec has integer values for seconds and nano seconds
    return (t->tv_sec + 1.0e-9 * (t->tv_nsec));
}
