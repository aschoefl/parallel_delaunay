#ifndef BUCKET_HPP
#define BUCKET_HPP


# include <mpi.h>
# include <random>
# include <vector>
# include <iterator>
# include <iostream>
# include <algorithm>
# include <memory>
# include <list>
# include <fstream>
# include <ctime>
# include "point.hpp"
# include "polygon.hpp"

using namespace std;


# define MAX_PNTS 10 // maximal points per bucket
# define MAX_PROC 500
# define TAG_NOTIFY 42
# define TAG_DATA 1
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
    static int P; // amount of processor R = P*P
    /* reveive buffer for MPI */
    static double buffer[MAX_PNTS*2+1];
    /* make sure there is only one root 
    achtung pfusch
    */
    static shared_ptr<Bucket> createRoot(int r){
        if (root == nullptr) {
            if (N%P) throw runtime_error("N must be divisible by P");
            if (P*P > MAX_PROC) 
                throw runtime_error("Maximal number of processors exceeded");
            int i = (r%P)*(N/P)+(N/P)/2;
            int j = (r/P)*(N/P)+(N/P)/2;
            cout << "create root with indices: ("<<i<<", "<<j<<")"<< endl;
            // srand(r*static_cast<unsigned int>(time(nullptr))); // set different seed for each processor
            srand(r*42); // set different seed for each processor
            Bucket* tmp = new Bucket(i,j);
            root = tmp->self;
            return tmp->self;
        } 
        cout << "root already exists" << endl;
        return root;
    };
    static shared_ptr<Bucket> root;
    static shared_ptr<Bucket> bb;
    static list<shared_ptr<Bucket>> buckets; // just for debugging

    shared_ptr<Bucket> operator() (int i, int j) const;

    void test(); //just for testing new functions

    inline void printNeighbours() const;
    /* don't just return false, because of downcasing*/ 
    inline bool isBnd( ) {return is_bnd;};
    void getIndex(int const dir, int& i, int& j);
    // void getIndex(Point const p, int& i, int& j);

    inline int dist(const shared_ptr<Bucket> other) const; 
    
    /*** get fuctions ***/
    int i() const { return ind_i;}
    int j() const { return ind_j;}
    int r() const {return (ind_j)/(N/P)*P + (ind_i)/(N/P);}
    // TODO: think about if return by value really is a good idea...
    vector<Point> getPoints();
    void printList();
    void doSomething();

/* TODO: make that nicer, not safe ? */
protected:
    Bucket() {}; // for boundary bucket

private:
    /*** private ctors and such ***/
    Bucket(int i, int j) : ind_i(i), ind_j(j) {
        self = move(shared_ptr<Bucket>(this));
        is_bnd = 0;
        addToList();
    };
    Bucket(const Bucket&); // deactivate copy-ctor
    Bucket& operator=(const Bucket&); // deactivate assign
    
    /*** private members ***/ 
    shared_ptr<Bucket> self; 
    /* global coordinates of points in order {cnt,x1,y1,..,xk,yk} 
        where cnt is the amount of points in the vector
    */
    vector<double> coordinates; 
    /* same points as point vector */
    vector<Point> points;

    int ind_i, ind_j; // global indices
    vector<shared_ptr<Bucket>> neighbours = {nullptr,nullptr,nullptr,nullptr,
        nullptr,nullptr,nullptr,nullptr }; 
    bool is_bnd = 1;
    bool running = 0;

    /*** private methods ***/
    void addPoint(double x, double y);
    void addBucket(int dir);
    void addToList();
    void fillCoordinates();
    void sendCoordinates(int destination, int no_bucket);
    void calculateDelauney();

};

class BoundaryBucket: private Bucket {
    inline bool isBnd( ) {return 1;};
};

#endif