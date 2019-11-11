
#ifndef MEDYAN_CAMKII_MathFunctions_h
#define MEDYAN_CAMKII_MathFunctions_h

#include <cmath>
#include <vector>

#include "common.h"

namespace mathfunc {

	/// Function to create an initial CaMKIIing point, given an
	/// initial normal vector and point.
	/// @param l - the distance of the camkii from the original point
	/// @return a vector describing the initial CaMKIIing position
	vector<double> camkiiProjection(const vector<double>& n, const vector<double>& p, double l);

}

#endif //MEDYAN_CAMKII_MathFunctions_h
