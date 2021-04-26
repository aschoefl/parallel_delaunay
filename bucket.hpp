#ifndef BUCKET_HPP
#define BUCKET_HPP


# include <mpi.h>
# include <vector>
# include <iterator>
# include <iostream>
# include <algorithm>
# include <memory>
# include <list>
#include <fstream>
# include "point.hpp"

using namespace std;

/*
NOTES
-----

*) using smart pointer because dtor was a problem,
   or I had memory leaks.

*) neighbours of bucket:
 3 | 2 | 1 
-----------
 4 |   | 0
-----------
 5 | 6 | 7

*/


class  Bucket {

public:
    static int N; // global amount of buckets is N*N
    static MPI_Status status; 
    static int P; // amount of processor R = P*P
    /* make sure there is only one root 
    achtung pfusch
    */
    static shared_ptr<Bucket> createRoot(int r){
        if (root == nullptr) {
            if (N%P) throw runtime_error("N must be divisible by P");
            int i = (r%P)*(N/P)+(N/P)/2;
            int j = (r/P)*(N/P)+(N/P)/2;
            cout << "create root with indices: ("<<i<<", "<<j<<")"<< endl;
            Bucket* tmp = new Bucket(i,j);
            root = tmp->self;
            return tmp->self;
        } 
        cout << "root already exists" << endl;
        return root;
    }
    int r() const {
        return (ind_j)/(N/P)*P + (ind_i)/(N/P);
    }
    static shared_ptr<Bucket> root;
    static shared_ptr<Bucket> bb;
    static list<Point> buckets; // just for debugging

    shared_ptr<Bucket> operator() (int i, int j) const;

    void test(); //just for testing new functions

    void setCoordinates(vector<double>&& vec) {
        coordinates = move(vec);
    };
    inline void printNeighbours() const;
    /* don't just return false, because of downcasing*/ 
    inline bool isBnd( ) {return is_bnd;};
    void getIndex(int dir, int& i, int& j);

    inline int dist(const shared_ptr<Bucket> other) const; 
    
    /*** get fuctions ***/
    int i() const { return ind_i;}
    int j() const { return ind_j;}
    // int N() const {return N;};

    /* returned by value, to be compatible with set function */
    vector<double> getCoordinates() const { return coordinates; };
    void printList();

/* TODO: make that nicer, not safe ? */
protected:
    Bucket() {}; // for boundary bucket

private:
    /*** private ctors and such ***/
    Bucket(int i, int j) : ind_i(i), ind_j(j) {
        self = move(shared_ptr<Bucket>(this));
        is_bnd = 0;
        addToList();
        // cout << "in ctor" << endl;
    };
    Bucket(const Bucket&); // deactivate copy-ctor
    Bucket& operator=(const Bucket&); // deactivate assign
    
    /*** private members ***/ 
    shared_ptr<Bucket> self; 
    /* global coordinates of points in order {x1,y1,..,xk,yk} */
    vector<double> coordinates; 
    int ind_i, ind_j; // global indices
    vector<shared_ptr<Bucket>> neighbours = {nullptr,nullptr,nullptr,nullptr,
        nullptr,nullptr,nullptr,nullptr }; 
    bool is_bnd = 1;
    vector<Point> points;

    /*** private methods ***/
    void addBucket(int dir);
    void addToList();
};

class BoundaryBucket: private Bucket {
    inline bool isBnd( ) {return 1;};
};

#endif