
#include "CaMKIIMathFunctions.h"

#include <cmath>
#include <vector>
#include <math.h>

#include "MathFunctions.h"
#include "Rand.h"

namespace mathfunc {

	vector<double> camkiiProjection(const vector<double> &n, const vector<double> &p, double l) {

		//get random permutation from p
		vector<double> r = {p[0] + Rand::randDouble(-1, 1),
							p[1] + Rand::randDouble(-1, 1),
							p[2] + Rand::randDouble(-1, 1)};

		//construct vector z which is r-p
		auto z = twoPointDirection(p, r);

		//construct u and v, which creates an orthogonal set n, u, v
		auto u = crossProduct(n, z);
		auto v = crossProduct(n, u);

		normalize(u);
		normalize(v);

		//find random point on circle defining the camkiiing point
		double thetaRandom = Rand::randDouble(0, 2 * M_PI);

		// the first point of CaMKII projection
		vector<double> cp;
		cp.push_back(p[0] + l * (u[0] * cos(thetaRandom) + v[0] * sin(thetaRandom)));
		cp.push_back(p[1] + l * (u[1] * cos(thetaRandom) + v[1] * sin(thetaRandom)));
		cp.push_back(p[2] + l * (u[2] * cos(thetaRandom) + v[2] * sin(thetaRandom)));

		return cp;
	}

}
