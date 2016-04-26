#include <cstdlib>
#include <iostream>
#include <cmath>
#include <omp.h>
#include "core.h"
#include "util.h"

#include <unistd.h>

#define L 200

#define PI 3.14159265359

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

int main()
{
	omp_set_num_threads(1);
	double *T = make_grid(L);
	double *b = make_grid(L);

	GRID_GET(b, L/4, L/4, L) = 1;
	//GRID_GET(b, 3*L/4, L/4, L) = 1;
	//GRID_GET(b, L/4, 3*L/4, L) = 1;
	GRID_GET(b, 3*L/4, 3*L/4, L) = 1;

	set_boundary(b, L, horizontal_boundary, vertical_boundary, horizontal_boundary, vertical_boundary);
	conjugate_gradient(T, b, L);
	dump_grid(T, L, "out.csv");
}

