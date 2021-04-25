/*
to execute: 
mpiexec -n 4  delauney 
*/ 
# include <iostream>
# include <cstdlib>
# include <mpi.h>
# include "polygon.hpp"
# include "point.hpp"
# include "bucket.hpp"
using namespace std;

// int main ( int argc, char *argv[] );


int main ( int argc, char *argv[] ) {

    /* local variables */
    int P, p, I; // no processes, rank, no local elements
    int finished = 0;
    MPI_Status status;

    // Polygon poly(0,0); 

    // poly.addPoint(-1,0);
    // poly.addPoint(1,0);
    // poly.addPoint(2,0);
    // poly.addPoint(0,-1);
    // poly.addPoint(0.5,-0.5);
    // poly.addPoint(1,-1);
    // poly.addPoint(0,1);

    // cout << poly << endl;

    // exit(0);

    /* Initialize MPI */
    // MPI_Init(&argc, &argv);
    // MPI_Comm_size(MPI_COMM_WORLD, &P);
    // MPI_Comm_rank(MPI_COMM_WORLD, &p);


    int N = 10; 

    shared_ptr<Bucket> root = Bucket::createRoot(N,5,5);
    root->test();

    exit(0);

    vector<PointBase> A[N][N];

    /* generate uniform example mesh */
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++) {
            A[i][j].push_back(PointBase((double)i/(double)N+1.0/(double)N/2.0, (double)j/(double)N+1.0/(double)N/2.0));
            // cout << "A["<<i<<","<<j<<"] = ";
            // for (auto tmp : A[i][j])
            //     cout << tmp;
            // cout << endl;
        }
    }

    // PointBase a(0,3);
    // PointBase b(1,3);
    // cout << (a+b)/2 << endl;

    /* choose one point in middle for testing
       for now test only the method for a larger polygon (in diagonal dir)
       TODO: implement spiral search as in paper
    */
    int i = N/2;
    int j = N/2;
    Polygon poly(A[i][j][0]);

    auto add = [&poly](const PointBase& p) {poly.addPoint(p);};
    vector<int> dir = {-1,1};
    for (int dir_i : dir) {
        for (int dir_j : dir) {
            auto incr = 1;
            bool stop = 0;
            while (!stop){

                // if (!(i+incr*dir_i==N || i+incr*dir_i==0 || j+incr*dir_j==N || i+incr*dir_j==0) ){
                //     if (dir_i==-1 && dir_j==-1) A[i+incr*dir_i][j+incr*dir_j].pop_back();
                // } // just for testing!

                /* add ghost point if out of bounds */
                if (i+incr*dir_i==N || i+incr*dir_i==0 || j+incr*dir_j==N || i+incr*dir_j==0) {
                    stop = 1;
                    add(PointBase(max(dir_i, 0), max(dir_j,0)));
                }
                /* add points (see sketch) and go in diagonal direction if cell is empty*/
                else if (!A[i+incr*dir_i][j+incr*dir_j].empty()) {
                    stop = 1;
                    for_each(A[i+incr*dir_i][j+incr*dir_j].begin(), A[i+incr*dir_i][j+incr*dir_j].end(), add);
                }
                else incr++;
            }
        }
    }
    poly.calculateVeroni();

    cout << A[i][j][0] << poly << endl;
    
    // MPI_Finalize();
    exit(0);
}


