# include "point.hpp"

bool Point::operator== (const Point& other) const {
    return x==other.x && y == other.y;
}

bool Point::operator!= (const Point& other) const {
    return !(*this == other);
}

bool Point::operator< (const Point& other) const {
    if (x == other.x) return y<other.y;
    return x<other.x;
}
bool Point::operator<= (const Point& other) const {
    if (*this == other) return 1;
    return *this < other;
}

bool Point::operator> (const Point& other) const {
    return !(*this <= other);
}

bool Point::operator>= (const Point& other) const {
    return !(*this < other);
}


std::ostream &operator<<(std::ostream &os, const Point &p) {
  os << "(" << p.x << ", " << p.y << ")";
  return os;
}