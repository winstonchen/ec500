#include <cstdlib>
#include <iostream>
#include <cmath>
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


/* === RED-BLACK ITERATION === */

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


/* === CONJUGATE GRADIENT METHOD === */

#define TOL 1e-4

static void apply_A(const double *T, double *T_out, const int L)
{
	const int max = L-1;

	for (int x = 1; x < max; x++) {
		for (int y = 1; y < max; y++) {
			GRID_GET(T_out, x, y, L) = 4*GRID_GET(T, x, y, L) -
			                           (GRID_GET(T, x+1, y, L) +
			                            GRID_GET(T, x-1, y, L) +
			                            GRID_GET(T, x, y+1, L) +
			                            GRID_GET(T, x, y-1, L));
		}
	}
}

static double dot_prod(const double *A, const double *B, const int L)
{
	double result = 0;
	for (int x = 0; x < L; x++) {
		for (int y = 0; y < L; y++) {
			result += GRID_GET(A, x, y, L)*GRID_GET(B, x, y, L);
		}
	}
	return result;
}

static double get_residual(const double *T, const double *b, double *r, const double alpha, const int L)
{
	apply_A(T, r, L);
	double res_mag = 0;
	for (int x = 0; x < L; x++) {
		for (int y = 0; y < L; y++) {
			const double res = GRID_GET(b, x, y, L) - alpha*GRID_GET(r, x, y, L);
			GRID_GET(r, x, y, L) = res;
			res_mag += res*res;
		}
	}
	return sqrt(res_mag);
}

void conjugate_gradient(double *T, const double *b, const int L)
{
	double *r = make_grid(L);
	double *r_new = make_grid(L);
	double *p = make_grid(L);
	double *A_p = make_grid(L);

	get_residual(T, b, r, 1, L);
	copy_grid(r, p, L);

	int k = 0;  // iteration number
	UNUSED(k);
	do {
		apply_A(p, A_p, L);
		const double alpha = dot_prod(r, r, L)/dot_prod(p, A_p, L);

		for (int x = 0; x < L; x++) {
			for (int y = 0; y < L; y++) {
				GRID_GET(T, x, y, L) += alpha*GRID_GET(p, x, y, L);
			}
		}

		const double res = get_residual(p, r, r_new, alpha, L);

		if (res < TOL)
			break;

		const double beta = dot_prod(r_new, r_new, L)/dot_prod(r, r, L);

		for (int x = 0; x < L; x++) {
			for (int y = 0; y < L; y++) {
				GRID_GET(p, x, y, L) = GRID_GET(r_new, x, y, L) + beta*GRID_GET(p, x, y, L);
			}
		}

		copy_grid(r_new, r, L);
		++k;
	} while (1);

	delete[] r;
	delete[] r_new;
	delete[] p;
	delete[] A_p;
}


