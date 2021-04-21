#ifndef BUCKET_HPP
#define BUCKET_HPP

# include <vector>
# include <iterator>
# include <iostream>
# include <algorithm>
# include <memory>
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
    /* make sure there is only one root*/
    static shared_ptr<Bucket> createRoot(int N_, int i, int j){
        Bucket* tmp = new Bucket(i,j);
        tmp->N = N_; 
        return tmp->self;
    }

    void setCoordinates(vector<double>&& vec) {
        coordinates = move(vec);
    };

    /* returned by value, to be compatible with set function */
    vector<double> getCoordinates(void) const {
        return coordinates;
    };

    void newBucket(int dir);
    int isCorner();

private: 
    Bucket();
    Bucket(int i, int j) : ind_i(i), ind_j(j) {
        self = move(shared_ptr<Bucket>(this));
        // cout << "in ctor" << endl;
    };
    Bucket(const Bucket&);

    shared_ptr<Bucket> self; 
    Bucket& operator=(const Bucket&);
    /* global coordinates of points in order {x1,y1,..,xk,yk} */
    vector<double> coordinates; 
    int ind_i, ind_j; // global indices
    void newToCorner(int diag);
    void addBucket(int dir);
    /* TODO: Make that nicer */
    vector<shared_ptr<Bucket>> neighbours = {nullptr,nullptr,nullptr,nullptr,
        nullptr,nullptr,nullptr,nullptr }; 
};

#endif