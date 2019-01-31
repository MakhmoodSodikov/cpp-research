#include <iostream>
#include <iomanip>
#include <vector>
#include <list>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <string>
#include <functional>
#include <unordered_map>
 
using namespace std;
 
 
template <typename Type>
Type abs (const Type &value) {
    return value < 0 ? -value : value;
}
 
 
// Vector 3D: ----------------------------------------------------------------------------------------------------
 
template <typename Type>
struct Vector3D {
    Type x, y, z;
 
    Vector3D (const Type &x = 0, const Type &y = 0, const Type &z = 0) {
        this->x = x;
        this->y = y;
        this->z = z;
    }
};
 
template <typename Type>
Vector3D<Type> operator + (const Vector3D<Type> &firstVector, const Vector3D<Type> &secondVector) {
    return Vector3D<Type>(firstVector.x + secondVector.x, firstVector.y + secondVector.y, firstVector.z + secondVector.z);
}
 
template <typename Type>
Vector3D<Type> operator * (const Type &number, const Vector3D<Type> &vec) {
    return Vector3D<Type>(number * vec.x, number * vec.y, number * vec.z);
}
 
template <typename Type>
Vector3D<Type> operator * (const Vector3D<Type> &vec, const Type &number) {
    return number * vec;
}
 
template <typename Type>
Vector3D<Type> operator / (const Type &number, const Vector3D<Type> &vec) {
    return vec / number;
}
 
template <typename Type>
Vector3D<Type> operator / (const Vector3D<Type> &vec, const Type &number) {
    return 1/number * vec;
}
 
template <typename Type>
Vector3D<Type> operator - (const Vector3D<Type> &firstVector, const Vector3D<Type> &secondVector) {
    return Vector3D<Type>(firstVector.x + -1 * secondVector.x, firstVector.y -1 * secondVector.y, firstVector.z + -1 * secondVector.z);
}
 
template <typename Type>
Type dotProduct (const Vector3D<Type> &firstVector, const Vector3D<Type> &secondVector) {
    return firstVector.x*secondVector.x + firstVector.y*secondVector.y + firstVector.z*secondVector.z;
}
 
template <typename Type>
Vector3D<Type> crossProduct (const Vector3D<Type> &firstVector, const Vector3D<Type> &secondVector) {
    return Vector3D<Type>(firstVector.y*secondVector.z - firstVector.z*secondVector.y, firstVector.x*secondVector.z - firstVector.z*secondVector.x, firstVector.x*secondVector.y - firstVector.y*secondVector.x);
}
 
template <typename Type>
Type getLength (const Vector3D<Type> &vec) {
    return sqrt(dotProduct(vec, vec));
}
 
 
//Point3D<Type>  3D: ----------------------------------------------------------------------------------------------------
 
template <typename Type>
struct Point3D {
    typedef Type CoordinateType;
    Type x, y, z;
 
    Point3D (const Type &x = 0, const Type &y = 0, const Type &z = 0) {
        this->x = x;
        this->y = y;
        this->z = z;
    }
};
 
template <typename Type>
bool operator == (const Point3D<Type> &firstPoint, const Point3D<Type> &secondPoint) {
    return (firstPoint.x == secondPoint.x && firstPoint.y == secondPoint.y && firstPoint.z == secondPoint.z);
}
 
template <typename Type>
bool operator != (const Point3D<Type> &firstPoint, const Point3D<Type> &secondPoint) {
    return !(firstPoint == secondPoint);
}
 
template <typename Type>
bool operator < (const Point3D<Type> &firstPoint, const Point3D<Type> &secondPoint) {
    return firstPoint.x < secondPoint.x || (firstPoint.x == secondPoint.x && (firstPoint.y < secondPoint.y || (firstPoint.z == secondPoint.z && firstPoint.z < secondPoint.z)));
}
 
template <typename Type>
bool operator <= (const Point3D<Type> &firstPoint, const Point3D<Type> &secondPoint) {
    return firstPoint < secondPoint || firstPoint == secondPoint;
}
 
template <typename Type>
bool operator > (const Point3D<Type> &firstPoint, const Point3D<Type> &secondPoint) {
    return !(firstPoint < secondPoint) && firstPoint != secondPoint;
}
 
template <typename Type>
bool operator >= (const Point3D<Type> &firstPoint, const Point3D<Type> &secondPoint) {
    return firstPoint > secondPoint || firstPoint == secondPoint;
}
 
template <typename Type>
Vector3D<Type> operator - (const Point3D<Type> &firstPoint, const Point3D<Type> &secondPoint) {
    return Vector3D<Type>(firstPoint.x - secondPoint.x, firstPoint.y - secondPoint.y, firstPoint.z - secondPoint.z);
}
 
template <typename Type>
Point3D<Type> operator + (const Point3D<Type> &firstPoint, const Point3D<Type> &secondPoint) {
    return Point3D<Type>(firstPoint.x + secondPoint.x, firstPoint.y + secondPoint.y, firstPoint.z + secondPoint.z);
}
 
template <typename Type>
Point3D<Type> operator + (const Point3D<Type> &point, const Vector3D<Type> &vec) {
    return Point3D<Type>(point.x + vec.x, point.y + vec.y, point.z + vec.z);
}
 
 
// Plane 3D: ----------------------------------------------------------------------------------------------------
 
template <class Type>
class Plane3D {
    private:
        Vector3D<Type> normal;
        Point3D<Type> point;
 
    public:
        Plane3D (const Point3D<Type> &firstPoint, const Point3D<Type> &secondPoint, const Point3D<Type> &thirdPoint) : point(firstPoint) {
            Vector3D<Type> firstVector = secondPoint - firstPoint;
            Vector3D<Type> secondVector = thirdPoint - firstPoint;
            normal = crossProduct(firstVector, secondVector);
            normal = normal / getLength(normal);
        }
 
        const Vector3D<Type> &getNormal () {
            return normal;
        }
};
 
template <typename Type>
ostream& operator << (ostream &os, const Point3D<Type> &point) {
    os << point.x << " " << point.y << " " << point.z;
    return os;
}
 
 
// Convex Hull 3D: ----------------------------------------------------------------------------------------------------
 
template<typename Type>
class ConvexHull3D {
    private:
        class Face {
            private:
                vector<Point3D<Type> > points;
       
            public:
                void addPoint (const Point3D<Type>  &point) {
                    points.push_back(point);
                }
 
                void swap () {
                    std::swap(points[1], points[2]);
                }
 
                const vector<Point3D<Type> > &getPoints () const {
                    return points;
                }
        };
 
        class Edge {
            private:
                pair<Point3D<Type>, Point3D<Type> > edge;
 
            public:
                struct EdgeHash {
                    size_t operator() (const Edge &edge) const {
                        return hash<Type>()(edge.edge.first.x) ^ hash<Type>()(edge.edge.first.y) ^ hash<Type>()(edge.edge.first.z) ^hash<Type>()(edge.edge.second.x) ^ hash<Type>()(edge.edge.second.y) ^ hash<Type>()(edge.edge.second.z);
                    }
                };
 
                Edge (const Point3D<Type> &firstPoint, const Point3D<Type> &secondPoint) {
                    edge = make_pair(min<Point3D<Type> >(firstPoint, secondPoint), max<Point3D<Type> >(firstPoint, secondPoint));
                }
 
                bool operator == (const Edge &otherEdge) const {
                    return edge.first == otherEdge.edge.first && edge.second == otherEdge.edge.second;
                }
        };
 
        list<Face> faces;
 
    public:
        typedef Type CoordinateType;
 
        inline bool isPointOnFace (const Point3D<Type> &point, const Face &face) {
            for (const Point3D<Type> &facePoint : face.getPoints())
                if (point == facePoint)
                    return true;
 
            return false;
        }
 
        char alph[26] = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'};
 
        void build (const vector<Point3D<Type> > &points) {
            // Добавляем первую грань:
            faces.resize(1);
 
            while (faces.begin()->getPoints().size() != 3) {
                size_t minPointIndex = numeric_limits<size_t>::max();
 
                for (size_t pointIndex = 0; pointIndex < points.size(); pointIndex ++) {
                    if (isPointOnFace(points[pointIndex], faces.back()))
                        continue;
 
                    if (minPointIndex == numeric_limits<size_t>::max() || points[pointIndex] < points[minPointIndex])
                        minPointIndex = pointIndex;
                }
 
                faces.begin()->addPoint(points[minPointIndex]);
            }
 
            for (size_t pointIndex = 0; pointIndex < points.size(); pointIndex ++) {
                if (!isPointOnFace(points[pointIndex], faces.back())) {
                    Plane3D<Type> plane(faces.back().getPoints()[0], faces.back().getPoints()[1], faces.back().getPoints()[2]);
                    Vector3D<Type> vec = faces.back().getPoints()[0] - points[pointIndex];
 
                    if (dotProduct(vec, plane.getNormal()) < 0)
                        prev(faces.end())->swap();
                   
                    break;
                }
            }
 
            unordered_map<Edge, unsigned short, typename Edge::EdgeHash> usedEdges;
            usedEdges[Edge(faces.back().getPoints()[0], faces.back().getPoints()[1])] = 1;
            usedEdges[Edge(faces.back().getPoints()[0], faces.back().getPoints()[2])] = 1;
            usedEdges[Edge(faces.back().getPoints()[1], faces.back().getPoints()[2])] = 1;
 
 
            // Заворачивание:
            while (true) {
                size_t maxPointIndex = numeric_limits<size_t>::max();
   
                for (size_t pointIndex = 0; pointIndex < points.size(); pointIndex ++) {
                    Edge firstEdge = Edge(points[pointIndex], *prev(faces.back().getPoints().end(), 1));
                    Edge secondEdge = Edge(points[pointIndex], *prev(faces.back().getPoints().end(), 2));
   
                    if (isPointOnFace(points[pointIndex], faces.back()) || usedEdges[firstEdge] >= 2 || usedEdges[secondEdge] >= 2)
                        continue;
 
                    if (maxPointIndex == numeric_limits<size_t>::max())
                        maxPointIndex = pointIndex;
 
                    else {
                        Plane3D<Type> prevPlane = Plane3D<Type>(faces.back().getPoints()[0], faces.back().getPoints()[1], faces.back().getPoints()[2]);
                        Plane3D<Type> maxPlane = Plane3D<Type>(points[maxPointIndex], *prev(faces.back().getPoints().end(), 1), *prev(faces.back().getPoints().end(), 2));
                        Plane3D<Type> newPlane = Plane3D<Type>(points[pointIndex], *prev(faces.back().getPoints().end(), 1), *prev(faces.back().getPoints().end(), 2));
                       
                        bool isNewPlaneBetterThanTheOldOne = getLength(crossProduct(newPlane.getNormal(), prevPlane.getNormal())) <= getLength(crossProduct(maxPlane.getNormal(), prevPlane.getNormal()));
   
                        if (isNewPlaneBetterThanTheOldOne)
                            maxPointIndex = pointIndex;
                    }
                }
 
                if (maxPointIndex == numeric_limits<size_t>::max())
                    break;
 
                Face newFace;
                newFace.addPoint(*prev(faces.back().getPoints().end(), 1));
                newFace.addPoint(*prev(faces.back().getPoints().end(), 2));
                newFace.addPoint(points[maxPointIndex]);
                faces.push_back(newFace);
 
                Edge firstEdge = Edge(faces.back().getPoints()[0], faces.back().getPoints()[1]);
                Edge secondEdge = Edge(faces.back().getPoints()[1], faces.back().getPoints()[2]);
                Edge thirdEdge = Edge(faces.back().getPoints()[2], faces.back().getPoints()[0]);
   
                usedEdges[firstEdge] ++;
                usedEdges[secondEdge] ++;
                usedEdges[thirdEdge] ++;
 
                cout << endl << endl << "face" << endl;
                cout << newFace.getPoints()[0] << endl << newFace.getPoints()[1] << endl << newFace.getPoints()[2] << endl;
            }
        }
};
 
 
int main () {
    size_t pointsCount;
    cin >> pointsCount;
 
    vector<Point3D<long double> > points(pointsCount);
 
    for (Point3D<long double> &point: points)
        cin >> point.x >> point.y >> point.z;
 
    ConvexHull3D<long double> convexHull;
    convexHull.build(points);
 

    return 0;
}