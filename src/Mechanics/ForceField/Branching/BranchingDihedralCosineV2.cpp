
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "BranchingDihedralCosineV2.h"
#include "BranchingDihedral.h"

#include "BranchingPoint.h"
#include "Bead.h"

#include "MathFunctions.h"
#include "Cylinder.h"

using namespace mathfunc;

//This version uses a the following vectors to determine dihedral angles.
// b1 = c4 - c3;
// b2 = c2 - mp;
// b3 = c3 - mp;
// n1 = b1 x b2;
// n2 = b3 x b2;

//Forces calculated are on points c5,c2,mp,c3. c5 = c4-c3+c2;

//STEP1: Transformation of forces on c5, c2, mp, and c3 TO forces on c2, mp, c3 and c4.
// Fc2 = -dE/dc2 + -dE/dc5* dc5/dc2 = -dE/dc2 + -dE/dc5 = Fc2 + Fc5
// F_mp = Fmp
// Fc3  = -dE/dc3 + -dE/dc5 * dc5/dc3 = Fc3 -Fc5
// Fc4  = -dE/dc5 * dc5/dc4 = Fc5

// STEP2: Transofmration of forces from c2, mp, c3, c4 TO c1, c2, c3, c4
// Fc1 = (1-p) * Fmp
// Fc2 = Fc2 + p * Fmp
// Fc3 = Fc3
// Fc4 = Fc4

floatingpoint BranchingDihedralCosineV2::energy(
		floatingpoint *coord, size_t nint,
		unsigned int *beadSet, floatingpoint *kdih, floatingpoint *pos){

	int n = BranchingDihedral<BranchingDihedralCosineV2>::n;

	floatingpoint *coord1, *coord2, *coord3, *coord4, n1n2, U_i;
	floatingpoint *mp = new floatingpoint[3];
	floatingpoint *n1 = new floatingpoint[3];
	floatingpoint *n2 = new floatingpoint[3];

	floatingpoint U = 0.0;

	for(int i = 0; i < nint; i += 1) {

		coord1 = &coord[3 * beadSet[n * i]];
		coord2 = &coord[3 * beadSet[n * i + 1]];
		coord3 = &coord[3 * beadSet[n * i + 2]];
		coord4 = &coord[3 * beadSet[n * i + 3]];

		midPointCoordinate(mp, coord1, coord2, pos[i]);

		// b1 = c4 - c3;
		// b2 = c2 - mp;
		// b3 = c3 - mp;
		// n1 = b1 x b2;
		// n2 = b3 x b2;
		vectorProduct(n1, coord3, coord4, mp, coord2);
		vectorProduct(n2, mp, coord3, mp, coord2);

		normalizeVector(n1);
		normalizeVector(n2);
		n1n2 = dotProduct(n1, n2);

		U_i = kdih[i] * ( 1 - n1n2 );

		if(fabs(U_i) == numeric_limits<floatingpoint>::infinity()
		   || U_i != U_i || U_i < -1.0) {

			//set culprit and return
			BranchingInteractions::_branchingCulprit = BranchingPoint::getBranchingPoints()[i];

			return -1;
		}

		U += U_i;
	}

//	cout<<"CosineV2     "<<U<<endl;
	delete [] mp;
	delete [] n1;
	delete [] n2;

	return U;
}

void BranchingDihedralCosineV2::forces(
		floatingpoint *coord, floatingpoint *f, size_t nint,
		unsigned int *beadSet, floatingpoint *kdih, floatingpoint *pos){

	int n = BranchingDihedral<BranchingDihedralCosineV2>::n;

	double *coord1, *coord2, *coord3, *coord4;
	floatingpoint *f1, *f2, *f3, *f4;

	coord1 = new double[3];
	coord2 = new double[3];
	coord3 = new double[3];
	coord4 = new double[3];

	double *mp = new double[3];
	double *n1 = new double[3];
	double *n2 = new double[3];
	double *zero = new double[3]; zero[0] = 0; zero[1] = 0; zero[2] = 0;

	//@{
	double Fc1[3], Fc2[3], Fc3[3], Fc4[3], Fc5[3], Fmp[3];
	double vb1x, vb1y, vb1z, vb2x, vb2y, vb2z, vb3x, vb3y, vb3z;
	double vb2xm, vb2ym, vb2zm;
	double ax, ay, az, bx, by, bz;
	double rasq, rbsq, rgsq, rg, rginv, ra2inv, rb2inv, rabinv;
	double c,s;
	double p, df1, ddf1;
	double fg, fga, gaa, gbb, hg, hgb;
	double dtfx, dtfy, dtfz, dtgx, dtgy, dtgz, dthx, dthy, dthz, df, sx2, sy2,
			sz2;

	for(int i = 0; i < nint; i += 1) {

		for(int j = 0; j < 3; j++){
			coord1[j] = coord[3 * beadSet[n * i]+ j];
			coord2[j] = coord[3 * beadSet[n * i + 1]+ j];
			coord3[j] = coord[3 * beadSet[n * i + 2]+ j];
			coord4[j] = coord[3 * beadSet[n * i + 3]+ j];
		}

		f1 = &f[3 * beadSet[n * i]];
		f2 = &f[3 * beadSet[n * i + 1]];
		f3 = &f[3 * beadSet[n * i + 2]];
		f4 = &f[3 * beadSet[n * i + 3]];

		// b1 = c4 - c3;
		// b2 = c2 - mp;
		// b3 = c3 - mp;
		// n1 = b1 x b2;
		// n2 = b3 x b2;

		//@ LAMMPS version test
		// 1st bond
		vb1x = coord4[0] - coord3[0];
		vb1y = coord4[1] - coord3[1];
		vb1z = coord4[2] - coord3[2];

		// 2nd bond
		vb2x = coord2[0] - mp[0];
		vb2y = coord2[1] - mp[1];
		vb2z = coord2[2] - mp[2];

		vb2xm = vb2x;
		vb2ym = vb2y;
		vb2zm = vb2z;

		// 3rd bond
		vb3x = coord3[0] - mp[0];
		vb3y = coord3[1] - mp[1];
		vb3z = coord3[2] - mp[2];

		// c,s calculation

		ax = vb1y*vb2zm - vb1z*vb2ym;
		ay = vb1z*vb2xm - vb1x*vb2zm;
		az = vb1x*vb2ym - vb1y*vb2xm;
		bx = vb3y*vb2zm - vb3z*vb2ym;
		by = vb3z*vb2xm - vb3x*vb2zm;
		bz = vb3x*vb2ym - vb3y*vb2xm;

		//|b1x-b2|
		rasq = ax*ax + ay*ay + az*az;
		//|b3x-b2|
		rbsq = bx*bx + by*by + bz*bz;
		//|b2|^2
		rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm;
		//|b2|
		rg = sqrt(rgsq);

		rginv = ra2inv = rb2inv = 0.0;
		if (rg > 0) rginv = 1.0/rg;//1/|-b2|
		if (rasq > 0) ra2inv = 1.0/rasq;//1/|b1x-b2|
		if (rbsq > 0) rb2inv = 1.0/rbsq;//1/|b3x-b2|
		rabinv = sqrt(ra2inv*rb2inv);//1/|b1x-b2||b3x-b2|

		c = (ax*bx + ay*by + az*bz)*rabinv;//(b1x-b2).(b3x-b2)/|b1x-b2||b3x-b2|
		s = rg*rabinv*(ax*vb3x + ay*vb3y + az*vb3z);//|b2|((b1x-b2).b3)/|b1x-b2||b3x-b2|=cos<n1,b3>/sin<b1,b2>

		if (c > 1.0) c = 1.0;
		if (c < -1.0) c = -1.0;

		ddf1 = c;
		df1 = s;

		p = 1.0 - c;

		fg = vb1x*vb2xm + vb1y*vb2ym + vb1z*vb2zm;//b1.-b2
		hg = vb3x*vb2xm + vb3y*vb2ym + vb3z*vb2zm;//b3.-b2
		fga = fg*ra2inv*rginv;//(b1.-b2)/(|b1x-b2||-b2|)
		hgb = hg*rb2inv*rginv;//(b3.-b2)/(|b3x-b2||-b2|)
		gaa = -ra2inv*rg;//-|-b2|/|b1x-b2|
		gbb = rb2inv*rg;//|-b2|/|b3x-b2|

		dtfx = gaa*ax;//-|-b2|(b1x-b2_x)/|b1x-b2|=-|-b2|(n1cap_x)
		dtfy = gaa*ay;//-|-b2|(b1x-b2_y)/|b1x-b2|=-|-b2|(n1cap_y)
		dtfz = gaa*az;//-|-b2|(b1x-b2_z)/|b1x-b2|=-|-b2|(n1cap_z)
		dtgx = fga*ax - hgb*bx;//-|-b2||n1cap_x|-(b3.-b2)(b3x-b2_x)/(|b3x-b2||-b2|)=-|-b2||n1cap_x|-(b3.-b2)n2cap_x/|-b2|
		dtgy = fga*ay - hgb*by;//-|-b2||n1cap_y|-(b3.-b2)(b3x-b2_y)/(|b3x-b2||-b2|)=-|-b2||n1cap_y|-(b3.-b2)n2cap_y/|-b2|
		dtgz = fga*az - hgb*bz;//-|-b2||n1cap_z|-(b3.-b2)(b3x-b2_z)/(|b3x-b2||-b2|)=-|-b2||n1cap_z|-(b3.-b2)n2cap_z/|-b2|
		dthx = gbb*bx;//|-b2|(b3x-b2_x)/|b3x-b2|=|-b2|(n2cap_x)
		dthy = gbb*by;//|-b2|(b3x-b2_y)/|b3x-b2|=|-b2|(n2cap_y)
		dthz = gbb*bz;//|-b2|(b3x-b2_z)/|b3x-b2|=|-b2|(n2cap_z)

		df = -kdih[i] * df1;//-Kd.s

		sx2 = df*dtgx;//-Kd cos<n1,b3>/sin<b1,b2>(-|-b2||n1cap_x|-(b3.-b2)n2cap_x/|-b2|)
		sy2 = df*dtgy;
		sz2 = df*dtgz;

		//Forces calculated are on points c5,c2,mp,c3. c5 = c4-c3+c2;

		//I II III IV
		//c5 c2 mp c3

		Fc2[0] = df*dtfx;
		Fc2[1] = df*dtfy;
		Fc2[2] = df*dtfz;

		Fc5[0] = sx2 - Fc2[0];
		Fc5[1] = sy2 - Fc2[1];
		Fc5[2] = sz2 - Fc2[2];

		Fc3[0] = df*dthx;
		Fc3[1] = df*dthy;
		Fc3[2] = df*dthz;

		Fmp[0] = -sx2 - Fc3[0];
		Fmp[1] = -sy2 - Fc3[1];
		Fmp[2] = -sz2 - Fc3[2];

		//STEP1: Transformation of forces on c5, c2, mp, and c3 TO forces on c2, mp, c3 and c4.
		// Fc2 = Fc2 + Fc5
		// F_mp = Fmp
		// Fc3  = Fc3 - Fc5
		// Fc4  = Fc5

		Fc2[0] += Fc5[0];
		Fc2[1] += Fc5[1];
		Fc2[2] += Fc5[2];

		Fc3[0] += -Fc5[0];
		Fc3[1] += -Fc5[1];
		Fc3[2] += -Fc5[2];

		Fc4[0] = Fc5[0];
		Fc4[1] = Fc5[1];
		Fc4[2] = Fc5[2];

		// STEP2: Transofmration of forces from c2, mp, c3, c4 TO c1, c2, c3, c4
		// Fc1 = (1-p) * Fmp
		// Fc2 = Fc2 + p * Fmp
		// Fc3 = Fc3
		// Fc4 = Fc4

		double alpha = pos[i];

		Fc1[0] = (1-alpha) * Fmp[0];
		Fc1[1] = (1-alpha) * Fmp[1];
		Fc1[2] = (1-alpha) * Fmp[2];

		Fc2[0] += (alpha) * Fmp[0];
		Fc2[1] += (alpha) * Fmp[1];
		Fc2[2] += (alpha) * Fmp[2];

		//Add to force vector
		f1[0] += Fc1[0];
		f1[1] += Fc1[1];
		f1[2] += Fc1[2];

		f2[0] += Fc2[0];
		f2[1] += Fc2[1];
		f2[2] += Fc2[2];

		f3[0] += Fc3[0];
		f3[1] += Fc3[1];
		f3[2] += Fc3[2];

		f4[0] += Fc4[0];
		f4[1] += Fc4[1];
		f4[2] += Fc4[2];
		//@
		#ifdef CHECKFORCES_INF_NAN
		if(checkNaN_INF<floatingpoint>(f1, 0, 2)||checkNaN_INF<floatingpoint>(f2,0,2)
		   ||checkNaN_INF<floatingpoint>(f3,0,2) ||checkNaN_INF<floatingpoint>(f4,0,2)){
			cout<<"Branching Dihedral Force becomes infinite. Printing data "<<endl;

			auto b = BranchingPoint::getBranchingPoints()[i];
			auto cyl1 = b->getFirstCylinder();
			auto cyl2 = b->getSecondCylinder();
			cout<<"Cylinder IDs "<<cyl1->getId()<<" "<<cyl2->getId()<<" with cIndex "
			    <<cyl1->getStableIndex()<<" "<<cyl2->getStableIndex()<<" and bIndex "
			    <<cyl1->getFirstBead()->getStableIndex()<<" "
			    <<cyl1->getSecondBead()->getStableIndex()<<" "
			    <<cyl2->getFirstBead()->getStableIndex()<<" "
			    <<cyl2->getSecondBead()->getStableIndex()<<endl;

			cout<<"Printing coords"<<endl;
			cout<<coord1[0]<<" "<<coord1[1]<<" "<<coord1[2]<<endl;
			cout<<coord2[0]<<" "<<coord2[1]<<" "<<coord2[2]<<endl;
			cout<<coord3[0]<<" "<<coord3[1]<<" "<<coord3[2]<<endl;
			cout<<coord4[0]<<" "<<coord4[1]<<" "<<coord4[2]<<endl;
			cout<<"Printing force"<<endl;
			cout<<f1[0]<<" "<<f1[1]<<" "<<f1[2]<<endl;
			cout<<f2[0]<<" "<<f2[1]<<" "<<f2[2]<<endl;
			cout<<f3[0]<<" "<<f3[1]<<" "<<f3[2]<<endl;
			cout<<f4[0]<<" "<<f4[1]<<" "<<f4[2]<<endl;
			cout<<"Printing binary Coords"<<endl;
			printvariablebinary(coord1,0,2);
			printvariablebinary(coord2,0,2);
			printvariablebinary(coord3,0,2);
			printvariablebinary(coord4,0,2);
			cout<<"Printing binary Force"<<endl;
			printvariablebinary(f1,0,2);
			printvariablebinary(f2,0,2);
			printvariablebinary(f3,0,2);
			printvariablebinary(f4,0,2);
			exit(EXIT_FAILURE);
		}
		#endif
	}

	delete [] coord1;
	delete [] coord2;
	delete [] coord3;
	delete [] coord4;

	delete [] mp;
	delete [] n1;
	delete [] n2;
	delete [] zero;

}