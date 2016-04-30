#include <cstdlib>
#include <iostream>
#include <cmath>
#include <omp.h>
#include <mpi.h>
#include "core.h"
#include "util.h"

#include <unistd.h>

#define L 200

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

struct nodestruct {
	// Processor id
	int rank;

	// Position in the processors grid
	int row;
	int col;

	// Neighbors
	double north;
	double south;
	double west;
	double east;
};

struct nodestruct mynode;

/* Derived data types used to shift matrix data */
MPI_Datatype MPI_InteriorPointsRow;
MPI_Datatype MPI_InteriorPointsCol;

int main()
{
	// MPI calls return
	int ierr;
	// My rank in the global communicator 
	int myid;


	// omp_set_num_threads(1);
	double *T = make_grid(L);
	double *b = make_grid(L);

	GRID_GET(b, L/4, L/4, L) = 1;
	//GRID_GET(b, 3*L/4, L/4, L) = 1;
	//GRID_GET(b, L/4, 3*L/4, L) = 1;
	GRID_GET(b, 3*L/4, 3*L/4, L) = 1;

    /*
	set_boundary(b, L, horizontal_boundary, vertical_boundary, horizontal_boundary, vertical_boundary);
	conjugate_gradient(T, b, L);
    */

	set_boundary(T, L, std::cos, std::sin, std::cos, std::sin);
	iterate_red_black(T, b, L, 10000, 0.5);
	dump_grid(T, L, "out.csv");
}

