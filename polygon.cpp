# include "polygon.hpp"

// add points in counter clockwise order
void Polygon::addPoint(double p0, double p1) {
    Point p(this, p0, p1);
    auto pos = std::find_if(points.begin(), points.end(), [p](auto s) {
        return s < p;
    });
    points.insert(pos, move(p));
}

void Polygon::addPoint(const PointBase& pin) {
    Point p(this, pin);
    auto pos = std::find_if(points.begin(), points.end(), [p](auto s) {
        return s < p;
    });
    points.insert(pos, move(p));
}

std::ostream &operator<<(std::ostream &os, const Polygon &poly) {
    os << "[ ";
    for (const auto& p: poly.points)
        os << p << " ";
    os << "]";
    return os;
}


bool Polygon::Point::operator> (const Polygon::Point& other) const {
    Point a = *this-poly->c;
    Point b = other-poly->c;

    // to have a beginning (at angle == 0 in unit circle)
    if (a.y>=0 && b.y<0) return 1;
    if (a.y<0 && b.y>=0) return 0;
    if (a.y==0 && b.y==0) {
        if (a.x>0 || b.x>0) return a.x<b.x;
        else return a.x>b.x;
    }

    // check order 
    auto cross_prod = a.x*b.y-a.y*b.x;
    if (cross_prod > 0) return 1;
    if (cross_prod < 0) return 0;

    // if points have same angle take nearer one first
    if (a.x*a.x+a.y*a.y < b.x*b.x+b.y*b.y) return 1;
    
    return 0;
}

bool Polygon::Point::operator>= (const Polygon::Point& other) const {
    if (*this>other || *this==other) return 1;
    return 0;
}

bool Polygon::Point::operator< (const Polygon::Point& other) const {
    if (!(*this>=other)) return 1;
    return 0;
}

bool Polygon::Point::operator<= (const Polygon::Point& other) const {
    if (!(*this>other)) return 1;
    return 0;
}
