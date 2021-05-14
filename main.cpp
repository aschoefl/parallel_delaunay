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
    Bucket::N = 4;
    int r, R;

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &R);
    MPI_Comm_rank(MPI_COMM_WORLD, &r);

    auto P = sqrt(R);
    if (P-int(P) == 0.) Bucket::P = int(P);
    else throw runtime_error("sqrt of amount of processors must be int");

    // vector<double> test;
    // for (int i = 0; i<5; i++)
    //     test.push_back(i);

    // test.insert(test.begin()+5, 42);

    // test.erase(test.begin()+3,test.end());
    // while (test.size() > 3) {
    //     test.pop_back();
    // }

    // for (int i = 0; i<test.size(); i++)
    //     cout << test[i] << " ";
    // cout << endl;

    // cout << -1%6 <<endl;

    shared_ptr<Bucket> root = Bucket::createRoot(r);
    root->doSomething();
    cout << root->r() << ": all done" << endl;
    // root->test();

    MPI_Finalize();
    exit(0);
}


// 3: poly: points: [ (1, 1), (0, 1), (1, 0), ]
// voronoi: [ (0.5, 1.375), (0.354167, 0.354167), (1.375, 0.5), ]
// radii: [ 0.625, 0.73657, 0.625, ]
