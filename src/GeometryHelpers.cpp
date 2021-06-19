#include "GeometryHelpers.h"

// some utility functions that can be shared
namespace ref::utils {

typedef ROOT::Math::XYPoint Point;

// check if a square in a ring
inline bool in_ring(const Point &pt, double sx, double sy, double rmin, double rmax, double phmin, double phmax)
{
    if (pt.r() > rmax || pt.r() < rmin) {
        return false;
    }

    // check four corners
    std::vector<Point> pts {
        Point(pt.x() - sx/2., pt.y() - sy/2.),
        Point(pt.x() - sx/2., pt.y() + sy/2.),
        Point(pt.x() + sx/2., pt.y() - sy/2.),
        Point(pt.x() + sx/2., pt.y() + sy/2.),
    };
    for (auto &p : pts) {
        if (p.r() > rmax || p.r() < rmin || p.phi() > phmax || p.phi() < phmin) {
            return false;
        }
    }
    return true;
}

// check if a square is overlapped with the others
inline bool overlap(const Point &pt, double sx, double sy, const std::vector<Point> &pts)
{
    for (auto &p : pts) {
        auto pn = p - pt;
        if ((std::abs(pn.x()) < (1. - 1e-6)*sx) && (std::abs(pn.y()) < (1. - 1e-6)*sy)) {
            return true;
        }
    }
    return false;
}

// a helper function to recursively fill square in a ring
void add_rectangle(Point p, std::vector<Point> &res, double sx, double sy, double rmin, double rmax, double phmin, double phmax)
{
    // outside of the ring or overlapping
    if (!in_ring(p, sx, sy, rmin, rmax, phmin, phmax) || overlap(p, sx, sy, res)) {
        return;
    }

    res.emplace_back(p);

    // check adjacent squares
    add_rectangle(Point(p.x() + sx, p.y()), res, sx, sy, rmin, rmax, phmin, phmax);
    add_rectangle(Point(p.x() - sx, p.y()), res, sx, sy, rmin, rmax, phmin, phmax);
    add_rectangle(Point(p.x(), p.y() + sy), res, sx, sy, rmin, rmax, phmin, phmax);
    add_rectangle(Point(p.x(), p.y() - sy), res, sx, sy, rmin, rmax, phmin, phmax);
}

// fill squares
std::vector<Point> fillRectangles(Point ref, double sx, double sy, double rmin, double rmax, double phmin, double phmax)
{
    // convert (0, 2pi) to (-pi, pi)
    if (phmax > M_PI) {
        phmin -= M_PI;
        phmax -= M_PI;
    }
    // start with a seed square and find one in the ring
    // move to center
    ref = ref - Point(int(ref.x()/sx)*sx, int(ref.y()/sy)*sy);

    auto find_seed = [] (const Point &ref, int n, double sx, double sy, double rmin, double rmax, double phmin, double phmax) {
        for (int ix = -n; ix < n; ++ix) {
            for (int iy = -n; iy < n; ++iy) {
                Point pt(ref.x() + ix*sx, ref.y() + iy*sy);
                if (in_ring(pt, sx, sy, rmin, rmax, phmin, phmax)) {
                    return pt;
                }
            }
        }
        return ref;
    };

    std::vector<Point> res;
    ref = find_seed(ref, int(rmax/sx) + 2, sx, sy, rmin, rmax, phmin, phmax);
    add_rectangle(ref, res, sx, sy, rmin, rmax, phmin, phmax);
    return res;
}


} // ref::utils
