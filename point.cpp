# include "point.hpp"

bool PointBase::operator== (const PointBase& other) const {
    return x==other.x && y == other.y;
}

bool PointBase::operator!= (const PointBase& other) const {
    return !(*this == other);
}

bool PointBase::operator< (const PointBase& other) const {
    if (x == other.x) return y<other.y;
    return x<other.x;
}
bool PointBase::operator<= (const PointBase& other) const {
    if (*this == other) return 1;
    return *this < other;
}

bool PointBase::operator> (const PointBase& other) const {
    return !(*this <= other);
}

bool PointBase::operator>= (const PointBase& other) const {
    return !(*this < other);
}


std::ostream &operator<<(std::ostream &os, const PointBase &p) {
  os << "(" << p.x << ", " << p.y << ")";
  return os;
}