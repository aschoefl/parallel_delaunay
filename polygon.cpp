# include "polygon.hpp"

template <typename T> 
void circumcenter(T& ret, const T& a, const T& b, const T& c) {
    auto alpha1 = c.x-a.x; 
    auto beta1 = a.y-c.y;
    auto gamma1 = alpha1*((a+c)/2).x+beta1*((a+c)/2).y;
    auto alpha2 = c.x-b.x; 
    auto beta2 = b.y-c.y;
    auto gamma2 = alpha2*((b+c)/2).x+beta2*((b+c)/2).y;
    auto det = alpha1*beta2-alpha2*beta1;

    if (det == 0) {
        throw ("lines are parallel in circumcenter");
    } else {
        ret.x = (gamma1*beta2-gamma2*beta1)/det;
        ret.y = (gamma2*alpha1-gamma1*alpha2)/det;
    }
}

// add points in counter clockwise order
void Polygon::addPoint(double p0, double p1) {
    PointPoly p(this, p0, p1);
    auto pos = std::find_if(points.begin(), points.end(), [p](auto s) {
        return s < p;
    });
    points.insert(pos, move(p));
}

void Polygon::addPoint(const Point& pin) {
    PointPoly p(this, pin);
    auto pos = std::find_if(points.begin(), points.end(), [p](auto s) {
        return s < p;
    });
    points.insert(pos, move(p));
}

void Polygon::calculateVeroni(void) { // ToDo: spelling
    if (points.size() < 2) return;
    PointPoly cc(this);

    for(auto it = points.begin(); it!=points.end()-1; it++){
        auto next = it+1;
        circumcenter(cc, *it, *next, c);
        veroni.push_back(cc); // cc is copied
    }
    circumcenter(cc, points.front(), points.back(), c);
    veroni.push_back(cc);
}

std::ostream &operator<<(std::ostream &os, const Polygon &poly) {
    os << "points: [ ";
    for (const auto& p: poly.points)
        os << p << " ";
    os << "]" << endl << "veroni: [ " ;
    for (const auto& p: poly.veroni)
        os << p << " ";
    os << "]" ;
    return os;
}

bool Polygon::PointPoly::operator> (const Polygon::PointPoly& other) const {
    PointPoly a = *this-poly->c;
    PointPoly b = other-poly->c;

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

bool Polygon::PointPoly::operator>= (const Polygon::PointPoly& other) const {
    if (*this>other || *this==other) return 1;
    return 0;
}

bool Polygon::PointPoly::operator< (const Polygon::PointPoly& other) const {
    if (!(*this>=other)) return 1;
    return 0;
}

bool Polygon::PointPoly::operator<= (const Polygon::PointPoly& other) const {
    if (!(*this>other)) return 1;
    return 0;
}
