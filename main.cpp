/*
to execute: 
mpiexec -n 4  delauney 
*/ 
# include <iostream>
# include <cstdlib>
# include <mpi.h>
# include "polygon.hpp"

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
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &p);


    int N = 10;
    double A[N][N][2];

    /* generate uniform example mesh */
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++) {
            A[i][j][0] = (double)i/(double)N+1.0/(double)N/2.0;
            A[i][j][1] = (double)j/(double)N+1.0/(double)N/2.0;
            // printf("(A[%d,%d] = (%.3f, %.3f) \n",i,j, A[i][j][0], A[i][j][1]);
        }
    }

    /* choose one point in middle for testing
       for now test only the method for a larger polygon
    */
    int i = N/2;
    int j = N/2;

    for (int eps=-1; eps<2; eps+=2) {
        for (int del=-1; del<2; del+=2) {
            
        }
    }







    
    MPI_Finalize();
    exit(0);
}


