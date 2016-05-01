#include <cstddef>
#include <cstring>
#include <iostream>
#include <fstream>
#include "util.h"

double *make_grid(const int L)
{
	const size_t size = L*L;
	double *array = new double[size];
	memset(array, 0, size*sizeof(double));
	return array;
}

void copy_grid(const double *from, double *to, const int L)
{
	const double max = L*L;
	for (int i = 0; i < max; i++) {
		to[i] = from[i];
	}
}

int dump_grid(const double *grid, const int L, const char *filename)
{
	std::ofstream out(filename);

	if (!out.is_open()) {
		return 1;
	}

	for (int x = 0; x < L; x++) {
		for (int y = 0; y < L; y++) {
			out << x << ',' << y << ',' << GRID_GET(grid, x, y, L) << std::endl;
		}
	}

	out.close();
	return 0;
}

