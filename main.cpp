# include <iostream>
# include <cstdlib>
# include <mpi.h>
# include <cmath>
# include "point.hpp"
# include "bucket.hpp"
# include <chrono>

using namespace std;

/* saved in repo https://github.com/aschoefl/parallel_delaunay */ 

int main ( int argc, char *argv[] ) {

    Bucket::N = 9*16*4; 
    int r, R;

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &R);
    MPI_Comm_rank(MPI_COMM_WORLD, &r);

    /* make sure P divides N */
    auto P = sqrt(R);
    if (P-int(P) == 0.) Bucket::P = int(P);
    else throw runtime_error("sqrt of amount of processors must be int");

    /* run and time the program */
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
