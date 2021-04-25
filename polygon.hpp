#ifndef POLYGON_HPP
#define POLYGON_HPP

# include <vector>
# include <iterator>
# include <iostream>
# include <algorithm>
# include "point.hpp"
using namespace std;

class Polygon {
public:
    Polygon(double c0, double c1): c(this, c0,c1) {};
    Polygon(Point p): c(this, p) {};

    void addPoint(double p0, double p1);
    void addPoint(const Point& pin);

    void calculateVeroni(void);

    friend std::ostream &operator<<(std::ostream &os, const Polygon &poly);

    class PointPoly: public Point {
    public:
        PointPoly(Polygon* pol) : poly(pol) {};
        PointPoly(Polygon* pol, const Point& p): poly(pol), Point(p) {};
        PointPoly(Polygon* pol, double p0, double p1): poly(pol), Point(p0, p1) {};
        PointPoly(const PointPoly& p): poly(p.poly), Point(p.x, p.y){};
        
        bool operator< (const PointPoly& other) const;
        bool operator> (const PointPoly& other) const;
        bool operator<= (const PointPoly& other) const;
        bool operator>= (const PointPoly& other) const;

    private: 
        Polygon* poly; 
    };

private: 
    vector<PointPoly> points;
    vector<PointPoly> veroni;
    PointPoly c;
};

#endif