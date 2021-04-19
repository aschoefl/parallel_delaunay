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
    Polygon(PointBase p): c(this, p) {};

    void addPoint(double p0, double p1);
    void addPoint(const PointBase& pin);

    friend std::ostream &operator<<(std::ostream &os, const Polygon &poly);

    class Point: public PointBase {
    public:
        Point() {};
        Point(Polygon* pol, const PointBase& p): poly(pol), PointBase(p) {};
        Point(Polygon* pol, double p0, double p1): poly(pol), PointBase(p0, p1) {};
        Point(const Point& p): poly(p.poly), PointBase(p.x, p.y){};
        
        /* TODO: template instead of double implementation */
        friend Point operator-(Point lhs, const Point& rhs){
            lhs -= rhs;
            return lhs;
        } 
        friend Point operator+(Point lhs, const Point& rhs){
            lhs += rhs;
            return lhs;
        } 

        bool operator< (const Point& other) const;
        bool operator> (const Point& other) const;
        bool operator<= (const Point& other) const;
        bool operator>= (const Point& other) const;


    private: 
        Polygon* poly; 
    };

private: 
    vector<Point> points;
    Point c;
};

#endif