#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <cmath>
// #include <omp.h>
#include <mpi.h>
#include <math.h>
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


int main(int argc, char *argv[])
{
	int nproc;
	int rank;
	MPI_Status Stat[4];
	MPI_Request request[4];

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// omp_set_num_threads(1);
	/*
	double *T = make_grid(L);
	double *b = make_grid(L);
	*/

	double *T = make_grid(L/nproc + 2, L);
	double *b = make_grid(L/nproc + 2, L);

	/*
	double** T = new double*[siz/nproc + 2];
  	for (int i = 0; i < siz; i++) {
    	T[i] = new double[siz];
  	}

  	double** b = new double*[siz/nproc + 2];
  	for (int i = 0; i < siz; i++) {
    	b[i] = new double[siz];
  	}
  	*/

	// GRID_GET(b, L/4, L/4, L) = 1;
	//GRID_GET(b, 3*L/4, L/4, L) = 1;
	//GRID_GET(b, L/4, 3*L/4, L) = 1;
	// GRID_GET(b, 3*L/4, 3*L/4, L) = 1;
	GRID_GET(b, (L/nproc + 2)/4, (L/nproc + 2)/4, (L/nproc + 2)) = 1;
	GRID_GET(b, 3*(L/nproc + 2)/4, 3*(L/nproc + 2)/4, (L/nproc + 2)) = 1;

    /*
	set_boundary(b, L, horizontal_boundary, vertical_boundary, horizontal_boundary, vertical_boundary);
	conjugate_gradient(T, b, L);
    */

	// set_boundary(T, L, std::cos, std::sin, std::cos, std::sin);
	/*
	void set_boundary(double *T, const int L,
	              double (*f_north)(double),
	              double (*f_east)(double),
	              double (*f_south)(double),
	              double (*f_west)(double))
	*/
	// const int max = L-1;

	for (int x = 0; x < L / nproc; x++) {
		// GRID_GET(T, x, 0, L)   = f_north(((double)x)/L);
		GRID_GET(T, x + 1, L - 1, L / nproc) = cos((double)x + rank * L/nproc / L);
		// GRID_GET(T, x, max, L) = f_south(((double)x)/L);
		GRID_GET(T, x + 1, L - 1, L / nproc) = cos((double)x * + rank * L/nproc / L);
	}

	for (int y = 0; y < L; y++) {
		// GRID_GET(T, 0, y, L)   = f_east(((double)y)/L);
		GRID_GET(T, 0, y + 1, L) = sin((double) y + rank * L/nproc / L);
		// GRID_GET(T, max, y, L) = f_west(((double)y)/L);
		GRID_GET(T, L - 1, y + 1, L) = sin((double) y + rank * L/nproc / L);
	}


	// iterate_red_black(T, b, L, 10000, 0.5);
	// void iterate_red_black(double *T, const double *b, const int L, const int iters, const double a)
	// const int max = L-1;

	double* ghostf_s = new double[L];
  	double* ghostb_s = new double[L];
 	double* ghostf_r = new double[L];
  	double* ghostb_r = new double[L];

  	double iters = 10000;
  	double a = 0.5;

	for (int i = 0; i < iters; i++) {
		/* red */
		// for (int x = 1; x < max; x++) {
		for (int x = 1; x < L / nproc + 2 - 1; x++) {
			// for (int y = 1; y < max; y++) {
			for (int y = 1; y < L - 1; y++) {
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
		// for (int x = 1; x < max; x++) {
		for (int x = 1; x < L / nproc + 2 - 1; x++) {
			// for (int y = 1; y < max; y++) {
			for (int y = 1; y < L - 1; y++) {
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

	  	for (int i = 0; i < siz; i++) {
	  		ghostb_s[i] = GRID_GET(T, 1, i, L / nproc + 2 - 1);
	  		ghostf_s[i] = GRID_GET(T, L / nproc, i, L / nproc + 2 - 1);
	  	}

		int rf = (rank + 1) % nproc;
		int rb = (rank - 1) % nproc;
		if (rb < 0) rb = nproc - 1;

		int tag1 = 1;
		int tag2 = 2;

		cout<<ghostb_r[0]<<" "<<ghostf_s[0]<<endl;

		MPI_Isend(&ghostf_s, L, MPI_DOUBLE, rf , tag1, MPI_COMM_WORLD, &request[0]);
		cout<<"Proc " <<rank<<" sending to  "<<rf<<" with tag "<<tag1<<endl;
		MPI_Isend(&ghostb_s, L, MPI_DOUBLE, rb , tag2, MPI_COMM_WORLD, &request[1]);
		MPI_Irecv(&ghostf_r, L, MPI_DOUBLE, rf , tag2, MPI_COMM_WORLD, &request[2]);
		MPI_Irecv(&ghostb_r, L, MPI_DOUBLE, rb , tag1, MPI_COMM_WORLD, &request[3]);
		cout<<"Proc " <<rank<<" recv from  "<<rb<<" with tag "<<tag1<<endl;

		MPI_Waitall(4, request, Stat);

		cout<<"Received"<<endl;
    	cout<<ghostb_r[0]<<" "<<ghostf_s[0]<<endl;

	 	for (int i = 0; i < L; i++) {
	 		cout << rank << "in for loop " <<ghostf_r[i] << " " << ghostb_r[i] << endl;
	 		GRID_GET(T, 0, i, L/nproc + 2) = ghostf_r[i];
	 		GRID_GET(T, L/nproc + 1, i, L/nproc + 2) = ghostb_r[i];
	 	}

		cout<<"Copied"<<endl;
	}

	dump_grid(T, L/nproc + 2, L, "out.csv");
}

