#include "GeometryHelpers.h"

// some utility functions that can be shared
namespace ref::utils {

typedef ROOT::Math::XYPoint Point;

// check if a square in a ring
inline bool in_ring(const Point &pt, double side, double rmin, double rmax, double phmin, double phmax)
{
    if (pt.r() > rmax || pt.r() < rmin) {
        return false;
    }

    // check four corners
    std::vector<Point> pts {
        Point(pt.x() - side/2., pt.y() - side/2.),
        Point(pt.x() - side/2., pt.y() + side/2.),
        Point(pt.x() + side/2., pt.y() - side/2.),
        Point(pt.x() + side/2., pt.y() + side/2.),
    };
    for (auto &p : pts) {
        if (p.r() > rmax || p.r() < rmin || p.phi() > phmax || p.phi() < phmin) {
            return false;
        }
    }
    return true;
}

// check if a square is overlapped with the others
inline bool overlap(const Point &pt, double side, const std::vector<Point> &pts)
{
    for (auto &p : pts) {
        auto pn = (p - pt)/side;
        if ((std::abs(pn.x()) < 1. - 1e-6) && (std::abs(pn.y()) < 1. - 1e-6)) {
            return true;
        }
    }
    return false;
}

// a helper function to recursively fill square in a ring
void add_square(Point p, std::vector<Point> &res, double lside, double rmin, double rmax,
                double phmin, double phmax)
{
    // outside of the ring or overlapping
    if (!in_ring(p, lside, rmin, rmax, phmin, phmax) || overlap(p, lside, res)) {
        return;
    }

    res.emplace_back(p);

    // check adjacent squares
    add_square(Point(p.x() + lside, p.y()), res, lside, rmin, rmax, phmin, phmax);
    add_square(Point(p.x() - lside, p.y()), res, lside, rmin, rmax, phmin, phmax);
    add_square(Point(p.x(), p.y() + lside), res, lside, rmin, rmax, phmin, phmax);
    add_square(Point(p.x(), p.y() - lside), res, lside, rmin, rmax, phmin, phmax);
}

// fill squares
std::vector<Point> fillSquares(Point ref, double lside, double rmin, double rmax, double phmin, double phmax)
{
    // start with a seed square and find one in the ring
    // move to center
    ref = ref - Point(int(ref.x()/lside)*lside, int(ref.y()/lside)*lside);

    auto find_seed = [] (const Point &ref, int n, double side, double rmin, double rmax, double phmin, double phmax) {
        for (int ix = -n; ix < n; ++ix) {
            for (int iy = -n; iy < n; ++iy) {
                Point pt(ref.x() + ix*side, ref.y() + iy*side);
                if (in_ring(pt, side, rmin, rmax, phmin, phmax)) {
                    return pt;
                }
            }
        }
        return ref;
    };

    std::vector<Point> res;
    ref = find_seed(ref, int(rmax/lside) + 2, lside, rmin, rmax, phmin, phmax);
    add_square(ref, res, lside, rmin, rmax, phmin, phmax);
    return res;
}


} // ref::utils
