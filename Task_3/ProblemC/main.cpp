#include <algorithm>
#include <iterator>
#include <iostream>
#include <vector>
#include <cmath>


typedef double sp_type;
const static double epsilon = 0.001f;
const static double pi = 3.141592653589793238462643f;


class point {
public:
    point() = default;
    point(sp_type X, sp_type Y);
    sp_type distance(point& p);

    sp_type X;
    sp_type Y;

    point operator+(point &x);
    point operator-(point &x);

    point& operator-() {
        Y = -Y;
        X = -X;
        return *this;
    }
};


sp_type point::distance(point &p) {
    return sqrt(pow((p.X - X),2) + pow((p.Y - Y),2));
}


point::point(sp_type X, sp_type Y) : X(X), Y(Y){}


point point::operator+(point &x) {
    return {this->X+x.X, this->Y+x.Y};
}


point point::operator-(point &x) {
    return {this->X-x.X, this->Y-x.Y};
}


istream& operator >> (istream &is, point& x) {
    is >> x.X;
    is >> x.Y;
    return is;
}


ostream& operator << (ostream &is, point& x) {
    is << x.X;
    is << x.Y;
    return is;
}


ostream& operator << (ostream &is, vector<point>& x) {
    for (auto e: x) {
        is << e;
    }
    return is;
}


class polygon {
public:
    vector<point> fig;
    void push_back(point p) { fig.push_back(p); };
    int size() { return fig.size(); };
    polygon() = default;
    explicit polygon(vector<point>& f) : fig(f) {};
};


class Vector2D {
public:
    point End{};
    double length;
    explicit Vector2D(point end) : End(end) { point O(0,0); length = O.distance(end); };
    Vector2D(point begin, point end);
    double dotProduct(Vector2D& other);
    Vector2D& operator+(Vector2D& vec) {
         End = End + vec.End;
         return *this;
    }
};


double Vector2D::dotProduct(Vector2D &other) {
    return this->End.X*other.End.X + this->End.Y*other.End.Y;
}


Vector2D::Vector2D(point begin, point end) {
    point O(0,0);

    End.X = end.X - begin.X;
    End.Y = end.Y - begin.Y;

    length = End.distance(O);
}


