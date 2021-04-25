#ifndef POINT_HPP
#define POINT_HPP

# include <vector>
# include <iterator>
# include <iostream>
# include <algorithm>
using namespace std;

class Point {
public:
    Point() {};
    Point(double p0, double p1): x(p0), y(p1) {};
    Point(const Point& p): x(p.x), y(p.y) {};

    double x,y; 
    friend std::ostream &operator<<(std::ostream &os, const Point &p);

    bool operator== (const Point& other) const;
    bool operator!= (const Point& other) const;
    /* lexicographic order */
    bool operator< (const Point& other) const;
    bool operator<= (const Point& other) const;
    bool operator> (const Point& other) const;
    bool operator>= (const Point& other) const;
    
    /* Arithmetic operators here because of templates*/
    template <typename T>
    Point& operator-= (const T& other){
        x -= other.x;
        y -= other.y;
        return *this;
    }
    template <typename T>
    Point& operator+= (const T& other){
        x += other.x;
        y += other.y;
        return *this;
    }
    template <typename T>
    Point& operator/= (const T f){
        x /= f;
        y /= f;
        return *this;
    }
    template <typename T>
    Point& operator*= (const T f){
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
