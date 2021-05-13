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

    // // test.erase(test.begin()+3,test.begin()+6);
    // while (test.size() > 3) {
    //     test.pop_back();
    // }

    // for (int i = 0; i<test.size(); i++)
    //     cout << test[i] << " ";
    // cout << endl;

    // cout << -1%6 <<endl;

    shared_ptr<Bucket> root = Bucket::createRoot(r);
    root->doSomething();
    // root->test();

    MPI_Finalize();
    exit(0);
}


// points: [ (0.875, 0.875) (0.375, 0.875) (0.375, 0.375) (0.875, 0.375) ]
// voronoi: [ (0.625, 0.875) (0.375, 0.625) (0.625, 0.375) (0.875, 0.625) ]
// radii: [ 0.25 0.25 0.25 0.25 ]


// it 0
// points: [ (0.875, 0.875) (0.625, 0.875) (0.375, 0.875) (0.375, 0.375) (0.875, 0.375) ]
// voronoi: [ (0.5, 0.75) (0, -0.125) (0.375, 0.625) (0.625, 0.375) (0.875, 0.625) ]
// radii: [ 0.176777 0.976281 0.25 0.25 0.25 0.25 ]

