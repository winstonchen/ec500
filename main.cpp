#include <cstdlib>
#include <iostream>
#include <cmath>
#include <omp.h>
#include "core.h"
#include "util.h"

#include <unistd.h>

#define L 200

int main()
{
	omp_set_num_threads(1);
	double *T = make_grid(L);
	double *b = make_grid(L);
	set_boundary(T, L, std::cos, std::sin, std::cos, std::sin);
	iterate_red_black(T, b, L, 10000, 0.5);
	dump_grid(T, L, "out.csv");
}

