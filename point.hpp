#ifndef POINT_HPP
#define POINT_HPP

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

    double x,y; 
    friend std::ostream &operator<<(std::ostream &os, const PointBase &p);

    bool operator== (const PointBase& other) const;
    bool operator!= (const PointBase& other) const;
    /* lexicographic order */
    bool operator< (const PointBase& other) const;
    bool operator<= (const PointBase& other) const;
    bool operator> (const PointBase& other) const;
    bool operator>= (const PointBase& other) const;
    
    /* Arithmetic operators here because of templates*/
    template <typename T>
    PointBase& operator-= (const T& other){
        x -= other.x;
        y -= other.y;
        return *this;
    }
    template <typename T>
    PointBase& operator+= (const T& other){
        x += other.x;
        y += other.y;
        return *this;
    }
    template <typename T>
    PointBase& operator/= (const T f){
        x /= f;
        y /= f;
        return *this;
    }
    template <typename T>
    PointBase& operator*= (const T f){
        x *= f;
        y *= f;
        return *this;
    }
    template <typename S, typename T>
    friend S operator-(S lhs, const T& rhs){
        lhs -= rhs;
        return lhs;
    } 
    template <typename S, typename T>
    friend S operator+(S lhs, const T& rhs){
        lhs += rhs;
        return lhs;
    } 
    template <typename S, typename T>
    friend S operator/(S lhs, const T f){
        lhs /= f;
        return lhs;
    } 
    template <typename S, typename T>
    friend S operator*(S lhs, const T f){
        lhs /= f;
        return lhs;
    } 
};

# endif
