#include <cstdlib>
#include <iostream>
#include <omp.h>
#include "util.h"
#include "core.h"

void set_boundary(double *T, const int L,
	              double (*f_north)(double),
	              double (*f_east)(double),
	              double (*f_south)(double),
	              double (*f_west)(double))
{
	const int max = L-1;

	for (int x = 0; x < L; x++) {
		GRID_GET(T, x, 0, L)   = f_north(((double)x)/L);
		GRID_GET(T, x, max, L) = f_south(((double)x)/L);
	}

	for (int y = 0; y < L; y++) {
		GRID_GET(T, 0, y, L)   = f_east(((double)y)/L);
		GRID_GET(T, max, y, L) = f_west(((double)y)/L);
	}
}

void iterate_red_black(double *T, const double *b, const int L, const int iters, const double a)
{
	const int max = L-1;

	for (int i = 0; i < iters; i++) {
		/* red */
		#pragma omp parallel for collapse(2)
		for (int x = 1; x < max; x++) {
			for (int y = 1; y < max; y++) {
				if ((x + y)%2 == 0)
					continue;

				GRID_GET(T, x, y, L) = (1-a)*GRID_GET(T, x, y, L) +
				                       a*(GRID_GET(T, x+1, y, L) +
				                          GRID_GET(T, x-1, y, L) +
				                          GRID_GET(T, x, y+1, L) +
				                          GRID_GET(T, x, y-1, L))/4 +
				                       a*GRID_GET(b, x, y, L);
			}
		}

		/* black */
		#pragma omp parallel for collapse(2)
		for (int x = 1; x < max; x++) {
			for (int y = 1; y < max; y++) {
				if ((x + y)%2 == 1)
					continue;

				GRID_GET(T, x, y, L) = (1-a)*GRID_GET(T, x, y, L) +
				                       a*(GRID_GET(T, x+1, y, L) +
				                          GRID_GET(T, x-1, y, L) +
				                          GRID_GET(T, x, y+1, L) +
				                          GRID_GET(T, x, y-1, L))/4 +
				                       a*GRID_GET(b, x, y, L);
			}
		}
	}
}

