// g++ -o triangle-lattice triangleLattice.cpp -lrt
#include <iostream>
#include <cmath>
#include <time.h>
#include <cstdio>

using namespace std;

double timeInSeconds(struct timespec *t){
    // a timespec has integer values for seconds and nano seconds
    return (t->tv_sec + 1.0e-9 * (t->tv_nsec));
}

void triangleLattice(double *x[], double *b[], int N)
{  
    int iters = 0;
    bool done = false;
    double temp[N][N];
    for(int x = 0; x < N; x++){
        for(int y = 0; y < N; y++){
            temp[x][y] = 0.0;
        }
    }
    while (!done){
        double res = 0.0;
        iters++;
        for(int i = 1; i < N-1; i++){
        	for(int j = 1; j < N-1; j++){
        		//temp[i] = (0.5 * x[i]) + (0.5/2) * (x[i+1] + x[i-1]);
        		temp[i][j] = (b[i][j]) + ((1.0/6.0) * (x[i-1][j+1] + x[i-1][j] + x[i][j+1] + x[i][j-1] + x[i+1][j] + x[i+1][j-1]));
        	}
        }
        //residue
        for(int i = 0; i < N; i++){ 
        	for(int j = 0; j < N; j++){
        		res += ( (x[i][j]  - temp[i][j]) * (x[i][j] - temp[i][j]) );
            	x[i][j] = temp[i][j];
        	}
        }
        //cout << iters << endl;
        res = sqrt(res);
        //cout << "residue: " << res << endl;
        if(res < 1e-6) {
            done = true;
            //cout << "residue: " << res << endl;
        }

    }
    cout << "iterations " << iters << endl;
}

int main()
{
    struct timespec time1, time2;
	int N = 50; 
	double *x[N], *b[N];
    for(int a = 0; a < N; a++){
        x[a] = new double[N];
        b[a] = new double[N];
    }
    for(int i = 0; i < N; i++){
    	for(int j = 0; j < N; j++){
    		x[i][j] = 0.0;   // phi
            b[i][j] = 0.0;   // res
    	}
    }
	b[N/2][N/2] = 1.0;
  	const int len = 5;
    //int sizes[len] = {10, 50, 100, 500, 1000};
    clock_gettime(CLOCK_MONOTONIC, &time1);
    triangleLattice(x,b,N);
    clock_gettime(CLOCK_MONOTONIC, &time2);
    double testTime = timeInSeconds(&time2)-timeInSeconds(&time1);
    printf("time: %lf secs \n", testTime);
  	// for(int i = 0; i < N; i++){
  	// 	for(int j = 0; j < N; j++){
  	// 		cout << x[i][j] << ",";
  	// 	}
  	// 	cout << endl;
  	// }
	return 0;
}