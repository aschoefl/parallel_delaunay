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
// # include "polygon.hpp"

using namespace std;


# define MAX_PNTS 9 // maximal points per bucket
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

*) Polygon in same header because of include problems

*/
class Bucket;

class Polygon {
public:

    friend class Bucket;
    Polygon(double c0, double c1): c(this, c0,c1) {};
    Polygon(Point p): c(this, p) {};
    Polygon() {};
    // Polygon(Polygon&& poly): c(move(poly.c)), points(move(poly.points)), 
    //     voronoi(move(poly.voronoi)), radii(move(poly.radii)), V(move(poly.V)){
    //         cout << "in move ctor" << endl;
    //     };
    Polygon& operator=(Polygon&& poly){
        // cout << "move assign " << endl;
        c = move(poly.c);
        points = move(poly.points); 
        voronoi = move(poly.voronoi);
        radii = move(poly.radii);
        V = move(poly.V);
        return *this;
    }

    void addPoint(double p0, double p1);
    void addPoint(const Point& pin);

    void calculateVoronoi(void);

    void printPoints(string no);
    friend std::ostream &operator<<(std::ostream &os, const Polygon &poly);

    class PointPoly: public Point {
    public:
        PointPoly() {};
        PointPoly(Polygon* pol) : poly(pol) {};
        PointPoly(Polygon* pol, const Point& p): poly(pol), Point(p) {};
        PointPoly(Polygon* pol, double p0, double p1): poly(pol), Point(p0, p1) {};
        PointPoly(const PointPoly& p): poly(p.poly), Point(p.x, p.y){};
        // PointPoly(PointPoly&& p): poly(move(p.poly)), Point(p.x, p.y){};
        PointPoly& operator=(PointPoly&& p){ 
            x = p.x; y =p.y;
            poly = move(p.poly);
            return *this;
        }


        bool operator< (const PointPoly& other) const;
        bool operator> (const PointPoly& other) const;
        bool operator<= (const PointPoly& other) const;
        bool operator>= (const PointPoly& other) const;

    private: 
        Polygon* poly; 
    };

private: 
    /* points in polygon */ 
    vector<PointPoly> points;
    /* same order as polygon pnts, order is not to be changed*/ 
    vector<Point> voronoi; 
    /* same order as polygon pnts, order is not to be changed*/ 
    vector<double> radii;
    /* Dealauney Candidates*/ 
    vector<Point> V;
    /* center of polygon */ 
    PointPoly c;
};

template <typename T> 
bool circumcenter(T& ret, const T& a, const T& b, const T& c);
extern template bool circumcenter<Point>(Point& ret, const Point& a, const Point& b, const Point& c);


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
    inline bool indexOutOfBnds(int ii, int jj) const ;
    inline int dist(const shared_ptr<Bucket> other) const; 
    
    /*** get fuctions ***/
    int i() const { return ind_i;}
    int j() const { return ind_j;}
    int r() const {return (ind_j)/(N/P)*P + (ind_i)/(N/P);}
    int getPoints(vector<Point>& pnts, int stat);
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

    /*** private methods ***/
    void addPoint(double x, double y);
    void addBucket(int dir);
    void addToList();
    int fillCoordinates(int stat);
    void sendCoordinates(int destination, int no_bucket);
    int calculateDelauney(int step);
    int initialize(int step);

    /* variables for Delauney calculation */
    Polygon poly;
    int init_dir_i, init_dir_j, init_incr ;
    int di, dj, it;
    int test_ind;
};

class BoundaryBucket: private Bucket {
    inline bool isBnd( ) {return 1;};
};


#endif