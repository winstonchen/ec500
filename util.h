#pragma once

#define UNUSED(x) ((void)(x))

#define GRID_GET(grid, x, y, L) ((grid)[(x) + (y)*(L)])

double *make_grid(const int L);
int dump_grid(const double *grid, const int L, const char *filename);

