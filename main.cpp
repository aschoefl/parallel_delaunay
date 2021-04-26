/*
to execute: 
mpiexec -n 4  delauney 
*/ 
# include <iostream>
# include <cstdlib>
# include <mpi.h>
# include <cmath>
# include "polygon.hpp"
# include "point.hpp"
# include "bucket.hpp"
using namespace std;

// int main ( int argc, char *argv[] );


int main ( int argc, char *argv[] ) {

    /* initiate variables */
    Bucket::N = 10;
    int r, R;

    /* Initialize MPI */
    // MPI_Init(&argc, &argv);
    // MPI_Comm_size(MPI_COMM_WORLD, &R);
    // MPI_Comm_rank(MPI_COMM_WORLD, &r);

    R = 25;
    r = 11;

    auto P = sqrt(R);
    if (P-int(P) == 0.) Bucket::P = int(P);
    else throw runtime_error("sqrt of amount of processors must be int");

    shared_ptr<Bucket> root = Bucket::createRoot(r);
    root->test();
    
    // MPI_Finalize();
    exit(0);
}


