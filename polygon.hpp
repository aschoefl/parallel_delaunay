# include <vector>
# include <iterator>
# include <iostream>
# include <algorithm>
using namespace std;

class Polygon {
public:
    Polygon(double c0, double c1): c(this, c0,c1) {};

    void addPoint(double p0, double p1);

    friend std::ostream &operator<<(std::ostream &os, const Polygon &poly);

    class Point {
    public:
        Point() {};
        Point(Polygon* pol, double p0, double p1): poly(pol), x(p0), y(p1) {};
        Point(const Point& p): poly(p.poly), x(p.x), y(p.y) {};

        friend std::ostream &operator<<(std::ostream &os, const Point &p);

        Point& operator-= (const Point& other);
        Point& operator+= (const Point& other);
        
        friend Point operator-(Point lhs, const Point& rhs){
            lhs -= rhs;
            return lhs;
        } // not both const ref for chained application

        friend Point operator+(Point lhs, const Point& rhs){
            lhs += rhs;
            return lhs;
        } // not both const ref for chained application


        bool operator== (const Point& other) const;
        bool operator< (const Point& other) const;
        bool operator> (const Point& other) const;
        bool operator<= (const Point& other) const;
        bool operator>= (const Point& other) const;

    protected:
        double x,y;

    private: 
        Polygon* poly; 
    };

private: 
    vector<Point> points;
    Point c;
};