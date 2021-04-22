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

    void test(); //just for testing new functions


    void setCoordinates(vector<double>&& vec) {
        coordinates = move(vec);
    };


    void newBucket(int dir);
    int isCorner();
    inline int dist(const shared_ptr<Bucket> other) const; 
    
    /*** get fuctions ***/
    int i() const { return ind_i;}
    int j() const { return ind_j;}
    /* returned by value, to be compatible with set function */
    vector<double> getCoordinates() const { return coordinates; };


private: 

    /*** private ctors and such ***/
    Bucket();
    Bucket(int i, int j) : ind_i(i), ind_j(j) {
        self = move(shared_ptr<Bucket>(this));
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

    /*** private methods ***/
    void addToCorner(int diag);
    void addBucket(int dir);
    shared_ptr<Bucket> searchCorner(vector<int>& to_go);
};

#endif