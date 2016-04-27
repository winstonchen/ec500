#include <iostream>
#include <cmath>


using namespace std;

void triangleLattice(double *x, double *b, int max_iters, int N)
{  
    int iters = 0;
    bool done = false;
    double temp[N][N];
    while (!done){
        double res = 0.0;
        iters++;
        for(int ii = 0; ii < N; ii++) temp[ii] = 0.0;
        for(int i = 1; i < N; i++){
        	for(int j = 1; j < N; j++){
        		//temp[i] = (0.5 * x[i]) + (0.5/2) * (x[i+1] + x[i-1]);
        		temp[i][j] = (x[i][j]) + ((1/6) * (x[i-1][j+1] + x[i-1][j] + x[i][j+1] + x[i][j-1] + x[i+1][j] + x[i+1][j-1]));
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
    return iters;
}

int main()
{
	int N = 10000;
	int max_iters = 5;
	double x[N][N], b[N][N];
    for(int i = 0; i < N; i++){
    	for(int j = 0; j < N; j++){
    		x[i][j] = 0.0;   // phi
            b[i][j] = 0.0;   // res
    	}
    }
	x[N/2][N/2] = 1.0;
  	const int len = 5;
    //int sizes[len] = {10, 50, 100, 500, 1000};

  	for(int i = 0; i < N; i++){
  		for(int j = 0; j < N; j++){
  			cout << x[ii][jj] << "," << endl;
  		}
  		cout << endl;
  	}
	return 0;
}