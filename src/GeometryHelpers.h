#pragma once
#include <vector>
#include "Math/Point2D.h"

// some utility functions that can be shared
namespace ref::utils {

typedef ROOT::Math::XYPoint Point;

// fill squares in a ring
std::vector<Point> fillSquares(Point ref, double lside, double rmin, double rmax,
                               double phmin = -M_PI, double phmax = M_PI);

} // ref::utils
