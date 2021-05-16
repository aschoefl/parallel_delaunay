/*
to execute: 
mpiexec -n 4  delauney 
*/ 
# include <iostream>
# include <cstdlib>
# include <mpi.h>
# include <cmath>
# include "point.hpp"
# include "bucket.hpp"
# include <chrono>

using namespace std;


int main ( int argc, char *argv[] ) {

    /* initiate variables */    
    Bucket::N = 20;//9*16*4; 
    int r, R;

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &R);
    MPI_Comm_rank(MPI_COMM_WORLD, &r);

    auto P = sqrt(R);
    if (P-int(P) == 0.) Bucket::P = int(P);
    else throw runtime_error("sqrt of amount of processors must be int");

    cout << "P = " << Bucket::P << ", N = " << Bucket::N << endl;

    shared_ptr<Bucket> root = Bucket::createRoot(r);
    auto start = chrono::high_resolution_clock::now();
    root->doSomething();
    auto stop = chrono::high_resolution_clock::now();

    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);

    cout << "processor " << root->r() << " done in " << duration.count() << " microseconds" << endl;

    MPI_Finalize();
    exit(0);
}


// 3: poly: points: [ (1, 1), (0, 1), (1, 0), ]
// voronoi: [ (0.5, 1.375), (0.354167, 0.354167), (1.375, 0.5), ]
// radii: [ 0.625, 0.73657, 0.625, ]
