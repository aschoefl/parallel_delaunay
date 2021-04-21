#ifndef BUCKET_HPP
#define BUCKET_HPP

# include <vector>
# include <iterator>
# include <iostream>
# include <algorithm>
# include "point.hpp"
using namespace std;

/*
neighbours of bucket:

 3 | 2 | 1 
-----------
 4 |   | 0
-----------
 5 | 6 | 7

*/


class  Bucket {
public:
    static int N; // global amount of buckets is N*N
    /* create only one instance outside of class to avoid holes*/
    static Bucket& createRoot(int N_, int i, int j){
        static Bucket root(i,j);
        N = N_; 
        return root;
    }

    void setCoordinates(vector<double>&& vec) {
        coordinates = move(vec);
    };

    /* returned by value, to be compatible with set function */
    vector<double> getCoordinates(void) {
        return coordinates;
    };

    void newBucket(int dir);

    /*TODO: DESTRUCTOR!*/ 

private: 
    Bucket() {};
    Bucket(int i, int j) : ind_i(i), ind_j(j) {};
    // Bucket(int i, int j, vector<double>&& vec) :
    //     ind_i(i), ind_j(j), coordinates(move(vec)) {};

    Bucket(const Bucket&);
    Bucket& operator=(const Bucket&);
    /* global coordinates of points in order {x1,y1,..,xk,yk} */
    vector<double> coordinates; 
    int ind_i, ind_j; // global indices

    /* TODO: Make that nicer!! */
    vector<Bucket*> neighbours = {nullptr,nullptr,nullptr,nullptr,
        nullptr,nullptr,nullptr,nullptr }; 

    inline int inc(int dir, int incr){
        if (incr < 0) throw("invalid increment in inc, must be postive");
        return (dir+incr)%8;
    }



};

#endif