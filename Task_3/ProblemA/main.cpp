#include <iostream>
#include <vector>
#include <cmath>


typedef double sp_type;
const sp_type E = 0.00000001;


class Point {
public:
    Point() = default;
    sp_type distance(Point& p);

    sp_type X;
    sp_type Y;
    sp_type Z;

    Point operator+(Point& x);
    Point operator-(Point &x);
    Point getDelta(Point R);
};


sp_type Point::distance(Point &p) {
    return sqrt(pow((p.X - X),2) + pow((p.Y - Y),2) + pow((p.Z - Z),2));
}


Point Point::getDelta(Point R) {
    sp_type delta_x = (R.X - X) / 3;
    sp_type delta_y = (R.Y - Y) / 3;
    sp_type delta_z = (R.Z - Z) / 3;
    return {delta_x, delta_y, delta_z};
}


Point Point::operator+(Point &x) {
    return {this->X+x.X, this->Y+x.Y, this->Z+x.Z};
}


Point Point::operator-(Point &x) {
    return {this->X-x.X, this->Y-x.Y, this->Z-x.Z};
}


class Segment {
public:
    Point Begin;
    Point End;
};


std::istream &operator>>(std::istream &is, Segment &s) {
    is >> s.Begin.X >> s.Begin.Y >> s.Begin.Z;
    is >> s.End.X >> s.End.Y >> s.End.Z;
    return is;
}


std::ostream& operator << (std::ostream& is, Segment& s) {
    is << s.Begin.X << s.Begin.Y << s.Begin.Z;
    is << s.End.X << s.End.Y << s.End.Z;
    return is;
}


class Solve {
public:
    sp_type getAns();
private:
    void readPoints();
    Segment m_sSegm1;
    Segment m_sSegm2;
    sp_type ternarSearchForFixedCoord(Point other_point);
    sp_type iterPoints();
};


sp_type Solve::getAns() {
    this->readPoints();
    return this->iterPoints();
}


void Solve::readPoints() {
    std::cin >> m_sSegm1;
    std::cin >> m_sSegm2;
}


sp_type Solve::ternarSearchForFixedCoord(Point other_point) {
    Point left = m_sSegm2.Begin;
    Point right = m_sSegm2.End;

    while(right.distance(left) > E) {
        Point delta = left.getDelta(right);
        Point m1 = left+delta;
        Point m2 = right-delta;
        if (m1.distance(other_point) > m2.distance(other_point)) {
            left = m1;
        } else {
            right = m2;
        }
    }

    return left.distance(other_point);
}


sp_type Solve::iterPoints() {
    Point left = m_sSegm1.Begin;
    Point right = m_sSegm1.End;
    while(right.distance(left) > E) {
        Point delta = left.getDelta(right);
        Point m1 = left+delta;
        Point m2 = right-delta;
        if (ternarSearchForFixedCoord(m1) > ternarSearchForFixedCoord(m2)) {
            left = m1;
        } else {
            right = m2;
        }
    }
    return ternarSearchForFixedCoord(left);
}


int main() {
    Solve sl{};
    printf("%.7f", sl.getAns());
    return 0;
}
