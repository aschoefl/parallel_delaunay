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
    static shared_ptr<Bucket> bb;

    shared_ptr<Bucket> operator() (int i, int j) const;

    void test(); //just for testing new functions


    void setCoordinates(vector<double>&& vec) {
        coordinates = move(vec);
    };
    inline void printNeighbours() const;
    /* don't just return false, because of downcasing*/ 
    inline bool isBnd( ) {return is_bnd;};
    void getIndex(int dir, int& i, int& j);

    void newBucket(int dir);
    int isCorner();
    inline int dist(const shared_ptr<Bucket> other) const; 
    
    /*** get fuctions ***/
    int i() const { return ind_i;}
    int j() const { return ind_j;}
    /* returned by value, to be compatible with set function */
    vector<double> getCoordinates() const { return coordinates; };

/* TODO: make that nicer, not safe*/
protected:
    Bucket() {}; // for boundary bucket

private: 
    /*** private ctors and such ***/
    Bucket(int i, int j) : ind_i(i), ind_j(j) {
        self = move(shared_ptr<Bucket>(this));
        is_bnd = 0;
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

    /*** private methods ***/
    void addToCorner(int diag);
    void addBucket(int dir);
    shared_ptr<Bucket> searchCorner(vector<int>& to_go);
};

class BoundaryBucket: private Bucket {
    inline bool isBnd( ) {return 1;};
};

#endif