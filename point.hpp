# include <vector>
# include <iterator>
# include <iostream>
# include <algorithm>
using namespace std;

class PointBase {
public:
    PointBase() {};
    PointBase(double p0, double p1): x(p0), y(p1) {};
    PointBase(const PointBase& p): x(p.x), y(p.y) {};

    friend std::ostream &operator<<(std::ostream &os, const PointBase &p);

    PointBase& operator-= (const PointBase& other);
    PointBase& operator+= (const PointBase& other);
    
    /* TODO: use templates instead of double implementation */
    friend PointBase operator-(PointBase lhs, const PointBase& rhs){
        lhs -= rhs;
        return lhs;
    } 
    friend PointBase operator+(PointBase lhs, const PointBase& rhs){
        lhs += rhs;
        return lhs;
    } 

    bool operator== (const PointBase& other) const;
    bool operator!= (const PointBase& other) const;

// protected:
    double x,y;

};
