#pragma once

void set_boundary(double *T, const int L,
	              double (*f_north)(double),
	              double (*f_east)(double),
	              double (*f_south)(double),
	              double (*f_west)(double));

void iterate_red_black(double *T, const double *b, const int L, const int iters, const double a);

