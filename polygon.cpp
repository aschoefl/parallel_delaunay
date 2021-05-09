// # include "polygon.hpp"



// // template <typename T> 
// // void circumcenter(T& ret, const T& a, const T& b, const T& c) {
// //     auto alpha1 = c.x-a.x; 
// //     auto beta1 = a.y-c.y;
// //     auto gamma1 = alpha1*((a+c)/2).x+beta1*((a+c)/2).y;
// //     auto alpha2 = c.x-b.x; 
// //     auto beta2 = b.y-c.y;
// //     auto gamma2 = alpha2*((b+c)/2).x+beta2*((b+c)/2).y;
// //     auto det = alpha1*beta2-alpha2*beta1;

// //     if (det == 0) {
// //         throw runtime_error("lines are parallel in circumcenter");
// //     } else {
// //         ret.x = (gamma1*beta2-gamma2*beta1)/det;
// //         ret.y = (gamma2*alpha1-gamma1*alpha2)/det;
// //     }
// // }

// template <typename T> 
// bool circumcenter(T& ret, const T& a, const T& b, const T& c) {
//     auto d = 2 * (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y));
//     if (d == 0)
//         return false; 
//         // throw runtime_error("lines are parallel in circumcenter");
//     ret.x = ((a.x * a.x + a.y * a.y) * (b.y - c.y) + (b.x * b.x + b.y * b.y) * (c.y - a.y) + (c.x * c.x + c.y * c.y) * (a.y - b.y)) / d;
//     ret.y = ((a.x * a.x + a.y * a.y) * (c.x - b.x) + (b.x * b.x + b.y * b.y) * (a.x - c.x) + (c.x * c.x + c.y * c.y) * (b.x - a.x)) / d;
//     return true;
// }
// template bool circumcenter<Point>(Point& ret, const Point& a, const Point& b, const Point& c);

// // add points in counter clockwise order
// void Polygon::addPoint(double p0, double p1) {
//     addPoint(PointPoly(this, p0, p1));
// }

// void Polygon::addPoint(const Point& pin) {
//     PointPoly p(this, pin);
//     // if (std::find(points.begin(), points.end(), p) != points.end())
//     //     return; 
//     bool exists = false;
//     auto pos = std::find_if(points.begin(), points.end(), [p,&exists](auto s) {
//         if (s==p) {
//             exists=true;
//         }
//         return s < p;
//     });
//     if (exists) return; // if element alread exists don't add
//     points.insert(pos, move(p));
// }

// void Polygon::calculateVoronoi(void) { // ToDo: spelling

//     voronoi.clear();
//     radii.clear();
//     if (points.size() < 2) return;
//     PointPoly cc(this); // to be circumcenter

//     for(auto it = points.begin(); it!=points.end()-1; it++){
//         auto next = it+1;
//         circumcenter(cc, *it, *next, c);
//         voronoi.push_back(cc); // cc is copied
//         radii.push_back(Point::dist(cc,c));
//     }
//     circumcenter(cc, points.front(), points.back(), c);
//     voronoi.push_back(cc);
//     radii.push_back(Point::dist(cc,c));
// }

// std::ostream &operator<<(std::ostream &os, const Polygon &poly) {
//     os << "points: [ ";
//     for (const auto& p: poly.points)
//         os << p << " ";
//     os << "]" << endl << "voronoi: [ " ;
//     for (const auto& p: poly.voronoi)
//         os << p << " ";
//     os << "]" << endl << "radii: [ " ;
//     for (const auto& p: poly.radii)
//         os << p << " ";
//     os << "]";
//     return os;
// }

// void Polygon::printPoints(string no) {

//     ofstream myfile;
    
//     myfile.open("polyPoints"+no+".txt");
//     if (myfile.is_open()) {
//         myfile << "[" << c << ",";

//         for (const auto& p: points)
//                 myfile << p << ", ";
//         myfile << "]" << endl;
//         myfile.close();
//     }
//     else cout << "Unable to open file";

//     myfile.open("voronoiPoints"+no+".txt");
//     if (myfile.is_open()) {
//         myfile << "[";
//         for (const auto& p: voronoi)
//                 myfile << p << ", ";
//         myfile << "]" << endl;
//         myfile.close();
//     }
//     else cout << "Unable to open file";
// }

// bool Polygon::PointPoly::operator> (const Polygon::PointPoly& other) const {
//     PointPoly a = *this-poly->c;
//     PointPoly b = other-poly->c;

//     // to have a beginning (at angle == 0 in unit circle)
//     if (a.y>=0 && b.y<0) return 1;
//     if (a.y<0 && b.y>=0) return 0;
//     if (a.y==0 && b.y==0) {
//         if (a.x>0 || b.x>0) return a.x<b.x;
//         else return a.x>b.x;
//     }

//     // check order 
//     auto cross_prod = a.x*b.y-a.y*b.x;
//     if (cross_prod > 0) return 1;
//     if (cross_prod < 0) return 0;

//     // if points have same angle take nearer one first
//     if (a.x*a.x+a.y*a.y < b.x*b.x+b.y*b.y) return 1;
    
//     return 0;
// }

// bool Polygon::PointPoly::operator>= (const Polygon::PointPoly& other) const {
//     if (*this>other || *this==other) return 1;
//     return 0;
// }

// bool Polygon::PointPoly::operator< (const Polygon::PointPoly& other) const {
//     if (!(*this>=other)) return 1;
//     return 0;
// }

// bool Polygon::PointPoly::operator<= (const Polygon::PointPoly& other) const {
//     if (!(*this>other)) return 1;
//     return 0;
// }
