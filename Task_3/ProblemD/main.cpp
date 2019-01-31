#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <set>


double const INFTY = 1e9;


struct Point2D {
    double x;
    double y;
    double z;

    int number;
    Point2D *Next, *Previous;

    explicit Point2D(double x = 0,
                     double y = 0,
                     double z = 0,
                     int id = -1)
            : x(x), y(y), z(z)
            , number(id)
            , Next(nullptr)
            , Previous(nullptr) {}

    bool make_event() {
        if (Previous->Next != this) {
            replace_event();
            return true;
        } else {
            preset_event();
            return false;
        }
    }

    void preset_event() const {
        Previous->Next = Next;

        Next->Previous = Previous;
    }

    void replace_event() {
        Previous->Next = this;

        Next->Previous = this;
    }

    friend bool operator<(const Point2D& p1, const Point2D& p2) {
        return p1.x < p2.x;
    }

    friend Point2D operator-(const Point2D& p1, const Point2D& p2) {
        return Point2D(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
    }

    friend double vec(const Point2D &p1, const Point2D &p2) {
        return p1.x * p2.y - p1.y * p2.x;
    }
};


struct Face {
    int points[3];

    explicit Face(int fst, int sec, int thrd){
        points[0] = fst;
        points[1] = sec;
        points[2] = thrd;
    }

    friend bool operator<(const Face& f1, const Face& f2) {
        if (f1.points[0] < f2.points[0]) {
            return true;
        } else if (f1.points[0] > f2.points[0]) {
            return false;
        } else {
            if (f1.points[1] < f2.points[1]) {
                return true;
            } else if (f1.points[1] > f2.points[1]) {
                return false;
            } else {
                return f1.points[2] < f2.points[2];
            }
        }
    }
};

typedef Point2D Event;

void preset_iters(const std::vector<Event*> *evs, size_t p1, size_t p2, Point2D *&left, Point2D *&right,
                  std::vector<double> &next_t);

void main_loop(Point2D *&u, Point2D *&v);

void initialize(const std::vector<Event *> *Events_vec, const Point2D *u, const Point2D *v, size_t pnt1, size_t pnt2,
                Point2D *&left, Point2D *&right, std::vector<double> &next_t);

void get_min(double current_t, const std::vector<double> &next_t, int &min_i, double &min_t);

void main_for(const std::vector<Point2D> &pts, size_t medium, const std::vector<Event *> &answer, Point2D *fst_pnt,
              Point2D *sec_pnt);

void check_nm(const std::vector<Point2D> &pts, size_t medium, Point2D *current, Point2D *&fst_pnt, Point2D *&sec_pnt);

void main_while(const std::vector<Event *> *Events_vec, std::vector<Event *> &answer, size_t pnt1, size_t pnt2,
                double cur_trp, Point2D *&fst_pnt, Point2D *&sec_pnt);

void pts_eq_curr(const Point2D *current, Point2D *&fst_pnt, Point2D *&sec_pnt);

void check_every_event(std::vector<Face> &result, const std::vector<Event *> &events);

double get_sign(const Point2D *a, const Point2D *b, const Point2D *c) {
    if (a == nullptr || b == nullptr || c == nullptr) {
        return INFTY;
    }
    return (b->x - a->x) * (c->y - b->y) - (b->y - a->y) * (c->x - b->x);
}

double set_time(const Point2D *a, const Point2D *b, const Point2D *c) {
    if (a == nullptr || b == nullptr || c == nullptr) {
        return INFTY;
    }
    return ((b->x - a->x) * (c->z - b->z) - (b->z - a->z) * (c->x - b->x)) / get_sign(a, b, c);
}

std::vector<Event*> make_hull(std::vector<Point2D> &pts, size_t lft, size_t rght) {
    if (rght - lft == 1) {
        return std::vector<Event*>();
    }

    size_t medium = (lft + rght) / 2;

    std::vector<Event*> Events_vec[2] = {
            make_hull(pts, lft, medium),
            make_hull(pts, medium, rght)
    };

    std::vector<Event*> answer;

    Point2D* fst_pnt = &pts[medium - 1];
    Point2D* sec_pnt = &pts[medium];

    main_loop(fst_pnt, sec_pnt);

    size_t pnt1 = 0, pnt2 = 0;
    double cur_trp = -INFTY;

    main_while(Events_vec, answer, pnt1, pnt2, cur_trp, fst_pnt, sec_pnt);

    fst_pnt->Next = sec_pnt;
    sec_pnt->Previous = fst_pnt;

    main_for(pts, medium, answer, fst_pnt, sec_pnt);

    return answer;
}

void main_while(const std::vector<Event *> *Events_vec, std::vector<Event *> &answer, size_t pnt1, size_t pnt2,
                double cur_trp, Point2D *&fst_pnt, Point2D *&sec_pnt) {
    while(true) {
        Point2D* left = nullptr;
        Point2D* right = nullptr;
        std::vector<double> nxtT(6, INFTY);

        initialize(Events_vec, fst_pnt, sec_pnt, pnt1, pnt2, left, right, nxtT);

        int min_i;
        double min_t;

        get_min(cur_trp, nxtT, min_i, min_t);

        if (min_i == -1 || min_t >= INFTY) {
            break;
        }

        if (min_i == 0) {
            if (left->x < fst_pnt->x) {
                answer.push_back(left);
            }
            left->make_event();
            pnt1++;
        }

        if (min_i == 1) {
            if (right->x > sec_pnt->x) {
                answer.push_back(right);
            }
            right->make_event();
            pnt2++;
        }

        if (min_i == 2) {
            answer.push_back(sec_pnt);
            sec_pnt = sec_pnt->Next;
        }

        if (min_i == 3) {
            sec_pnt = sec_pnt->Previous;
            answer.push_back(sec_pnt);
        }

        if (min_i == 4) {
            answer.push_back(fst_pnt);
            fst_pnt = fst_pnt->Previous;
        }

        if (min_i == 5) {
            fst_pnt = fst_pnt->Next;
            answer.push_back(fst_pnt);
        }
        cur_trp = min_t;
    }
}

void main_for(const std::vector<Point2D> &pts, size_t medium, const std::vector<Event *> &answer, Point2D *fst_pnt,
              Point2D *sec_pnt) {
    for (int i = static_cast<int>(answer.size() - 1); i >= 0; i--) {
        Point2D* current = answer[i];
        if (current->x > fst_pnt->x && current->x < sec_pnt->x) {
            check_nm(pts, medium, current, fst_pnt, sec_pnt);
        } else {
            current->make_event();
            pts_eq_curr(current, fst_pnt, sec_pnt);
        }
    }
}

void pts_eq_curr(const Point2D *current, Point2D *&fst_pnt, Point2D *&sec_pnt) {
    if (current == fst_pnt) {
        fst_pnt = fst_pnt->Previous;
    }
    if (current == sec_pnt) {
        sec_pnt = sec_pnt->Next;
    }
}

void check_nm(const std::vector<Point2D> &pts, size_t medium, Point2D *current, Point2D *&fst_pnt, Point2D *&sec_pnt) {
    fst_pnt->Next = sec_pnt->Previous = current;
    current->Previous = fst_pnt;
    current->Next = sec_pnt;
    if (current->x <= pts[medium - 1].x) {
        fst_pnt = current;
    } else {
        sec_pnt = current;
    }
}

void get_min(double current_t, const std::vector<double> &next_t, int &min_i, double &min_t) {
    min_i= -1;
    min_t= INFTY;
    for (int i = 0; i < 6; i++) {
        if (next_t[i] > current_t && next_t[i] < min_t) {
            min_t = next_t[i];
            min_i = i;
        }
    }
}

void initialize(const std::vector<Event *> *Events_vec, const Point2D *u, const Point2D *v, size_t pnt1, size_t pnt2,
                Point2D *&left, Point2D *&right, std::vector<double> &next_t) {
    preset_iters(Events_vec, pnt1, pnt2, left, right, next_t);
    next_t[2] = set_time(u, v, v->Next);
    next_t[3] = set_time(u, v->Previous, v);
    next_t[4] = set_time(u->Previous, u, v);
    next_t[5] = set_time(u, u->Next, v);
}

void main_loop(Point2D *&u, Point2D *&v) {
    while(true) {
        if (get_sign(u, v, v->Next) < 0) {
            v = v->Next;
        } else if (get_sign(u->Previous, u, v) < 0) {
            u = u->Previous;
        } else {
            break;
        }
    }
}

void preset_iters(const std::vector<Event *> *evs, size_t p1, size_t p2, Point2D *&left, Point2D *&right,
                  std::vector<double> &next_t) {
    if (p1 < evs[0].size()) {
        left = evs[0][p1];
        next_t[0] = set_time(left->Previous, left, left->Next);
    }
    if (p2 < evs[1].size()) {
        right = evs[1][p2];
        next_t[1] = set_time(right->Previous, right, right->Next);
    }
}


std::vector<Face> build_hll(std::vector<Point2D> pnts) {
    size_t n = pnts.size();
    std::vector<Face> res;
    std::sort(pnts.begin(), pnts.end());
    std::vector<Event*> events = make_hull(pnts, 0, n);

    check_every_event(res, events);

    for (Point2D& p : pnts) {
        p.Next = nullptr;
        p.Previous = nullptr;
        p.z = -p.z;
    }

    events = make_hull(pnts, 0, n);
    for (Event* event : events) {
        Face current(event->Previous->number, event->number, event->Next->number);

        if (event->make_event()) {
            std::swap(current.points[0], current.points[1]);
        }

        res.push_back(current);
    }
    return res;
}

void check_every_event(std::vector<Face> &result, const std::vector<Event *> &events) {
    for (Event* event : events) {
        Face current(event->Previous->number, event->number, event->Next->number);
        if (!event->make_event()) {
            std::swap(current.points[0], current.points[1]);
        }
        result.push_back(current);
    }
}


bool accept(const Point2D &a, const Point2D &b, const Point2D &c) {
    return vec(b - a, c - b) > 0;
}


typedef std::pair<int, int> Edge;


void check_all_sites(const std::vector<Point2D> &sites, std::vector<Point2D> &hull) {
    for (const Point2D& site : sites) {
        while (hull.size() >= 2) {
            if (accept(hull[hull.size() - 2], hull.back(), site)) {
                break;
            }
            hull.pop_back();
        }
        hull.push_back(site);
    }
}

void check_pts_in_sides(const std::vector<Point2D> &sites, const std::vector<int> &edges_number,
                        const std::vector<bool> &is_in_planar_hull, int &in_pnts, int &sum_deg) {
    for (int i = 0; i < sites.size(); i++) {
        if (!is_in_planar_hull[i]) {
            sum_deg += edges_number[i];
            in_pnts++;
        }
    }
}

double count_num(std::vector<Point2D> &sites) {
    std::vector<Face> dim_h = build_hll(sites);
    std::set<Edge> edges;
    std::vector<int> edges_number(sites.size());
    std::vector<bool> is_in_planar_hull(sites.size(), false);

    for (Face& face : dim_h) {
        for (int j = 0; j < 3; j++) {
            Edge edge(face.points[j], face.points[(j + 1) % 3]);
            if (edge.first > edge.second) {
                std::swap(edge.first, edge.second);
            }
            edges.insert(edge);
        }
    }

    for (const Edge& edge : edges) {
        edges_number[edge.first]++;
        edges_number[edge.second]++;
    }

    std::sort(sites.begin(), sites.end());
    std::vector<Point2D> hull;

    check_all_sites(sites, hull);

    for (int i = sites.size() - 2, bottom = hull.size(); i >= 0; --i) {
        while (static_cast<int>(hull.size()) > bottom) {
            if (accept(hull[hull.size() - 2], hull.back(), sites[i])) {
                break;
            }
            hull.pop_back();
        }
        hull.push_back(sites[i]);
    }

    for (Point2D& i : hull) {
        is_in_planar_hull[i.number] = true;
    }

    int in_pnts = 0;

    int sum_deg = 0;

    check_pts_in_sides(sites, edges_number, is_in_planar_hull, in_pnts, sum_deg);

    if (in_pnts == 0) {
        return 0.0;
    } else {
        return static_cast<double>(sum_deg) / in_pnts;
    }
}


class Solve {
public:
    Solve(){
        freopen("input.txt", "r", stdin);
        std::vector<Point2D> pts;
        double x, y;

        get_pts(pts, x, y);

        std::cout << count_num(pts) << "\n";
    }

    void get_pts(std::vector<Point2D> &points, double& x, double& y) const {
        for (int i = 0; std::cin >> y >> x; i++) {
            Point2D p(x, y, x * x + y * y, i);
            points.push_back(p);
        }
    }
};


int main() {
    Solve sl;
    return 0;
}
