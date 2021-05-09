// #ifndef POLYGON_HPP
// #define POLYGON_HPP

// # include <vector>
// # include <iterator>
// # include <iostream>
// # include <algorithm>
// # include <fstream>
// # include <cmath>
// # include "point.hpp"
// # include "bucket.hpp"
// using namespace std;

// class Polygon {
// public:

//     friend class Bucket;
//     Polygon(double c0, double c1): c(this, c0,c1) {};
//     Polygon(Point p): c(this, p) {};
//     // Polygon(Polygon&& poly): c(move(poly.c)), points(move(poly.points)), 
//     //     voronoi(move(poly.voronoi)), radii(move(poly.radii)), V(move(poly.V)){
//     //         cout << "in move ctor" << endl;
//     //     };
//     Polygon& operator=(Polygon&& poly){
//         cout << "move assign " << endl;
//         c = move(poly.c);
//         points = move(poly.points); 
//         voronoi = move(poly.voronoi);
//         radii = move(poly.radii);
//         V = move(poly.V);
//         return *this;
//     }

//     void addPoint(double p0, double p1);
//     void addPoint(const Point& pin);

//     void calculateVoronoi(void);

//     void printPoints(string no);
//     friend std::ostream &operator<<(std::ostream &os, const Polygon &poly);

//     class PointPoly: public Point {
//     public:
//         PointPoly(Polygon* pol) : poly(pol) {};
//         PointPoly(Polygon* pol, const Point& p): poly(pol), Point(p) {};
//         PointPoly(Polygon* pol, double p0, double p1): poly(pol), Point(p0, p1) {};
//         PointPoly(const PointPoly& p): poly(p.poly), Point(p.x, p.y){};
//         // PointPoly(PointPoly&& p): poly(move(p.poly)), Point(p.x, p.y){};
//         PointPoly& operator=(PointPoly&& p){ 
//             x = p.x; y =p.y;
//             poly = move(p.poly);
//             return *this;
//         }


//         bool operator< (const PointPoly& other) const;
//         bool operator> (const PointPoly& other) const;
//         bool operator<= (const PointPoly& other) const;
//         bool operator>= (const PointPoly& other) const;

//     private: 
//         Polygon* poly; 
//     };

// private: 
//     /* points in polygon */ 
//     vector<PointPoly> points;
//     /* same order as polygon pnts, order is not to be changed*/ 
//     vector<Point> voronoi; 
//     /* same order as polygon pnts, order is not to be changed*/ 
//     vector<double> radii;
//     /* Dealauney Candidates*/ 
//     vector<Point> V;
//     /* center of polygon */ 
//     PointPoly c;
// };

// template <typename T> 
// bool circumcenter(T& ret, const T& a, const T& b, const T& c);
// extern template bool circumcenter<Point>(Point& ret, const Point& a, const Point& b, const Point& c);

// #endif