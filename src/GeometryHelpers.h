#pragma once
#include <vector>
#include "Math/Point2D.h"

// some utility functions that can be shared
namespace ref::utils {

typedef ROOT::Math::XYPoint Point;

// fill rectangles in a ring
std::vector<Point> fillRectangles(Point ref, double sx, double sy, double rmin, double rmax,
                                  double phmin = -M_PI, double phmax = M_PI);
// fill squares in a ring
inline std::vector<Point> fillSquares(Point ref, double size, double rmin, double rmax,
                                      double phmin = -M_PI, double phmax = M_PI)
{
    return fillRectangles(ref, size, size, rmin, rmax, phmin, phmax);
}

} // ref::utils
