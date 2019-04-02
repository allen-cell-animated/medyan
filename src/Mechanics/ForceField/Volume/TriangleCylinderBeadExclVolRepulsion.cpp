
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2017-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "TriangleCylinderBeadExclVolRepulsion.h"

#include "Triangle.h"
#include "Vertex.h"
#include "Bead.h"

#include "MathFunctions.h"
#include "SysParams.h"

using namespace mathfunc;

/**
 * Implements the volume exclusion force field of triangle and bead.
 * 
 * Please refer to the document files for the math derivation of this force
 * field. Developers of this force field are responsible of keeping the
 * related documents up to date.
 */

// Note: These functions require that the area of the triangle has already been calculated

double TriangleCylinderBeadExclVolRepulsion::energy(
    Vec3 c0, Vec3 c1, Vec3 c2, Vec3 cb, double area,
    double kExVol
) {

    //check if in same plane
    auto cp = cross(c1 - c0, c2 - c0);
    if(areEqual(dot(cp, cb - c0), 0.0)) {
        
        // slightly move point. Using negative number here to move "into the cell".
        cb -= cp * (0.01 / magnitude(cp));
    }
    
    double A = dot(c1 - c0, c1 - c0);
    double B = dot(c2 - c1, c2 - c1);
    double C = dot(c0 - cb, c0 - cb);
    double D = 2 * dot(c1 - c0, c0 - cb);
    double E = 2 * dot(c2 - c1, c0 - cb);
    double F = 2 * dot(c1 - c0, c2 - c1);

    double A1 = 2 * A*E - D*F;
    double A2 = 2 * B*D - 2 * A*E + (D - E)*F;
    double A3 = -4 * A*B - 2 * B*D + F*(E + F);

    double B1 = 4 * A*C - D*D;
    double B2 = 4 * A*C - D*D + 4 * B*C - E*E + 4 * C*F - 2 * D*E;
    double B3 = 4 * B*A - F*F + 4 * B*C - E*E + 4 * B*D - 2 * E*F;
    double BB1 = sqrt(B1);
    double BB2 = sqrt(B2);
	double BB3 = sqrt(B3);

    double C1 = 2 * A + D;
    double C2 = 2 * A + D + E + 2 * (B + F);
    double C3 = 2 * B + E + F;
    double D1 = D;
    double D2 = D + E;
    double D3 = E + F;

    double E1 = atan(C1 / BB1);
    double E2 = atan(C2 / BB2);
    double E3 = atan(C3 / BB3);
    double F1 = atan(D1 / BB1);
    double F2 = atan(D2 / BB2);
    double F3 = atan(D3 / BB3);

    double G1 = A1 / BB1;
    double G2 = A2 / BB2;
    double G3 = A3 / BB3;

    double numerator = G1*(E1 - F1) + G2*(E2 - F2) + G3*(E3 - F3);
    double denominator = B*D*D + A*(-4 * B*C + E*E) + F*(-D*E + C*F);

    double H = numerator / denominator;
    
    double energy = kExVol * 2 * area * H;
    
    return energy;
}

void TriangleCylinderBeadExclVolRepulsion::forces(
    double* f0, double* f1, double* f2, double* fb,
    Vec3 c0, Vec3 c1, Vec3 c2, Vec3 cb,
    double area, const Vec3& dArea0, const Vec3& dArea1, const Vec3& dArea2,
    double kExVol
) {
    
    //check if in same plane
    auto cp = cross(c1 - c0, c2 - c0);
    if(areEqual(dot(cp, cb - c0), 0.0)) {
        
        // slightly move point. Using negative number here to move "into the cell".
        cb -= cp * (0.01 / magnitude(cp));
    }
    
    double A = dot(c1 - c0, c1 - c0);
    double B = dot(c2 - c1, c2 - c1);
    double C = dot(c0 - cb, c0 - cb);
    double D = 2 * dot(c1 - c0, c0 - cb);
    double E = 2 * dot(c2 - c1, c0 - cb);
    double F = 2 * dot(c1 - c0, c2 - c1);

    // Derivative index is 0, 1, 2, b
    array<Vec3, 4> dA = {
        2 * (c0[0] - c1[0]), 2 * (c0[1] - c1[1]), 2 * (c0[2] - c1[2]),
        2 * (c1[0] - c0[0]), 2 * (c1[1] - c0[1]), 2 * (c1[2] - c0[2])
    };
    array<Vec3, 4> dB = {
        0, 0, 0,
        2 * (c1[0] - c2[0]), 2 * (c1[1] - c2[1]), 2 * (c1[2] - c2[2]),
        2 * (c2[0] - c1[0]), 2 * (c2[1] - c1[1]), 2 * (c2[2] - c1[2])
    };
    array<Vec3, 4> dC = {
        2 * (c0[0] - cb[0]), 2 * (c0[1] - cb[1]), 2 * (c0[2] - cb[2]),
        0, 0, 0,
        0, 0, 0,
        2 * (cb[0] - c0[0]), 2 * (cb[1] - c0[1]), 2 * (cb[2] - c0[2])
    };
    array<Vec3, 4> dD = {
        2 * (c1[0] + cb[0] - 2*c0[0]), 2 * (c1[1] + cb[1] - 2*c0[1]), 2 * (c1[2] + cb[2] - 2*c0[2]),
        2 * (c0[0] - cb[0]), 2 * (c0[1] - cb[1]), 2 * (c0[2] - cb[2]),
        0, 0, 0,
        2 * (c0[0] - c1[0]), 2 * (c0[1] - c1[1]), 2 * (c0[2] - c1[2])
    };
    array<Vec3, 4> dE = {
        2 * (c2[0] - c1[0]), 2 * (c2[1] - c1[1]), 2 * (c2[2] - c1[2]),
        2 * (cb[0] - c0[0]), 2 * (cb[1] - c0[1]), 2 * (cb[2] - c0[2]),
        2 * (c0[0] - cb[0]), 2 * (c0[1] - cb[1]), 2 * (c0[2] - cb[2]),
        2 * (c1[0] - c2[0]), 2 * (c1[1] - c2[1]), 2 * (c1[2] - c2[2])
    };
    array<Vec3, 4> dF = {
        2 * (c1[0] - c2[0]), 2 * (c1[1] - c2[1]), 2 * (c1[2] - c2[2]),
        2 * (c0[0] + c2[0] - 2*c1[0]), 2 * (c0[1] + c2[1] - 2*c1[1]), 2 * (c0[2] + c2[2] - 2*c1[2]),
        2 * (c1[0] - c0[0]), 2 * (c1[1] - c0[1]), 2 * (c1[2] - c0[2])
    };

    double A1 = 2 * A*E - D*F;
    double A2 = 2 * B*D - 2 * A*E + (D - E)*F;
    double A3 = -4 * A*B - 2 * B*D + F*(E + F);

    double B1 = 4 * A*C - D*D;
    double B2 = 4 * A*C - D*D + 4 * B*C - E*E + 4 * C*F - 2 * D*E;
    double B3 = 4 * B*A - F*F + 4 * B*C - E*E + 4 * B*D - 2 * E*F;
    double BB1 = sqrt(B1);
    double BB2 = sqrt(B2);
	double BB3 = sqrt(B3);

    double C1 = 2 * A + D;
    double C2 = 2 * A + D + E + 2 * (B + F);
    double C3 = 2 * B + E + F;
    double D1 = D;
    double D2 = D + E;
    double D3 = E + F;

    double E1 = atan(C1 / BB1);
    double E2 = atan(C2 / BB2);
    double E3 = atan(C3 / BB3);
    double F1 = atan(D1 / BB1);
    double F2 = atan(D2 / BB2);
    double F3 = atan(D3 / BB3);

    double G1 = A1 / BB1;
    double G2 = A2 / BB2;
    double G3 = A3 / BB3;

    double numerator = G1*(E1 - F1) + G2*(E2 - F2) + G3*(E3 - F3);
    double denominator = B*D*D + A*(-4 * B*C + E*E) + F*(-D*E + C*F);

    double H = numerator / denominator;
    
    array<Vec3, 4> dH = {};
    for(size_t dIdx = 0; dIdx < 4; ++dIdx) {
        for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
            double dAT = dA[dIdx][coordIdx];
            double dBT = dB[dIdx][coordIdx];
            double dCT = dC[dIdx][coordIdx];
            double dDT = dD[dIdx][coordIdx];
            double dET = dE[dIdx][coordIdx];
            double dFT = dF[dIdx][coordIdx];

            double dA1 = 2 * (dAT*E + A*dET) - (dDT*F + D*dFT);
            double dA2 = 2 * (dBT*D + B*dDT) - 2 * (dAT*E + A*dET) + ((dDT - dET)*F + (D - E)*dFT);
            double dA3 = -4 * (dAT*B + A*dBT) - 2 * (dBT*D + B*dDT) + (dFT*(E + F) + F*(dET + dFT));

            double dB1 = 4 * (dAT*C + A*dCT) - 2 * D*dDT;
            double dB2 = 4 * (dAT*C + A*dCT) - 2 * D*dDT + 4 * (dBT*C + B*dCT) - 2 * E*dET + 4 * (dCT*F + C*dFT) - 2 * (dDT*E + D*dET);
            double dB3 = 4 * (dBT*A + B*dAT) - 2 * F*dFT + 4 * (dBT*C + B*dCT) - 2 * E*dET + 4 * (dBT*D + B*dDT) - 2 * (dET*F + E*dFT);
            double dBB1 = dB1 / 2 / BB1;
            double dBB2 = dB2 / 2 / BB2;
            double dBB3 = dB3 / 2 / BB3;

            double dC1 = 2 * dAT + dDT;
            double dC2 = 2 * dAT + dDT + dET + 2 * (dBT + dFT);
            double dC3 = 2 * dBT + dET + dFT;
            double dD1 = dDT;
            double dD2 = dDT + dET;
            double dD3 = dET + dFT;

            double dE1 = (BB1*dC1 - C1*dBB1) / (B1 + C1*C1);
            double dE2 = (BB2*dC2 - C2*dBB2) / (B2 + C2*C2);
            double dE3 = (BB3*dC3 - C3*dBB3) / (B3 + C3*C3);
            double dF1 = (BB1*dD1 - D1*dBB1) / (B1 + D1*D1);
            double dF2 = (BB2*dD2 - D2*dBB2) / (B2 + D2*D2);
            double dF3 = (BB3*dD3 - D3*dBB3) / (B3 + D3*D3);

            double dG1 = (BB1*dA1 - A1*dBB1) / B1;
            double dG2 = (BB2*dA2 - A2*dBB2) / B2;
            double dG3 = (BB3*dA3 - A3*dBB3) / B3;

            double dNumerator = dG1*(E1 - F1) + G1*(dE1 - dF1) + dG2*(E2 - F2) + G2*(dE2 - dF2) + dG3*(E3 - F3) + G3*(dE3 - dF3);
            double dDenominator = dBT * D*D + B * 2 * D*dDT + dAT*(-4 * B*C + E*E) + A*(-4 * (dBT*C + B*dCT) + 2 * E*dET) + dFT*(-D*E + C*F) + F*(-(dDT*E + D*dET) + dCT*F + C*dFT);

            dH[dIdx][coordIdx] = (denominator*dNumerator - numerator*dDenominator) / (denominator*denominator);
        }
    }

    for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
        fb[coordIdx] -= kExVol * dH[3][coordIdx] * 2 * area;
        f0[coordIdx] -= kExVol * 2 * (dArea0[coordIdx] * H + dH[0][coordIdx] * area);
        f1[coordIdx] -= kExVol * 2 * (dArea1[coordIdx] * H + dH[1][coordIdx] * area);
        f2[coordIdx] -= kExVol * 2 * (dArea2[coordIdx] * H + dH[2][coordIdx] * area);
    }
}

Vec3 TriangleCylinderBeadExclVolRepulsion::loadForces(
    const Vec3& c0, const Vec3& c1, const Vec3& c2, const Vec3& coord,
    double area, double kExVol
) {
    Vec3 cb = coord;

    //check if in same plane
    auto cp = cross(c1 - c0, c2 - c0);
    if(areEqual(dot(cp, cb - c0), 0.0)) {
        
        // slightly move point. Using negative number here to move "into the cell".
        cb -= cp * (0.01 / magnitude(cp));
    }
    
    double A = dot(c1 - c0, c1 - c0);
    double B = dot(c2 - c1, c2 - c1);
    double C = dot(c0 - cb, c0 - cb);
    double D = 2 * dot(c1 - c0, c0 - cb);
    double E = 2 * dot(c2 - c1, c0 - cb);
    double F = 2 * dot(c1 - c0, c2 - c1);

    // Only consider derivative on the coordinate of the bead (cb)
    // Vec3 dA = {{}};
    // Vec3 dB = {{}};
    Vec3 dC = { 2 * (cb[0] - c0[0]), 2 * (cb[1] - c0[1]), 2 * (cb[2] - c0[2]) };
    Vec3 dD = { 2 * (c0[0] - c1[0]), 2 * (c0[1] - c1[1]), 2 * (c0[2] - c1[2]) };
    Vec3 dE = { 2 * (c1[0] - c2[0]), 2 * (c1[1] - c2[1]), 2 * (c1[2] - c2[2]) };
    Vec3 dF = { 2 * (c1[0] - c0[0]), 2 * (c1[1] - c0[1]), 2 * (c1[2] - c0[2]) };

    double A1 = 2 * A*E - D*F;
    double A2 = 2 * B*D - 2 * A*E + (D - E)*F;
    double A3 = -4 * A*B - 2 * B*D + F*(E + F);

    double B1 = 4 * A*C - D*D;
    double B2 = 4 * A*C - D*D + 4 * B*C - E*E + 4 * C*F - 2 * D*E;
    double B3 = 4 * B*A - F*F + 4 * B*C - E*E + 4 * B*D - 2 * E*F;
    double BB1 = sqrt(B1);
    double BB2 = sqrt(B2);
	double BB3 = sqrt(B3);

    double C1 = 2 * A + D;
    double C2 = 2 * A + D + E + 2 * (B + F);
    double C3 = 2 * B + E + F;
    double D1 = D;
    double D2 = D + E;
    double D3 = E + F;

    double E1 = atan(C1 / BB1);
    double E2 = atan(C2 / BB2);
    double E3 = atan(C3 / BB3);
    double F1 = atan(D1 / BB1);
    double F2 = atan(D2 / BB2);
    double F3 = atan(D3 / BB3);

    double G1 = A1 / BB1;
    double G2 = A2 / BB2;
    double G3 = A3 / BB3;

    double numerator = G1*(E1 - F1) + G2*(E2 - F2) + G3*(E3 - F3);
    double denominator = B*D*D + A*(-4 * B*C + E*E) + F*(-D*E + C*F);

    // Only consider derivative on the coordinate of the bead (cb)
    Vec3 dH {};
    for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
        double dCT = dC[coordIdx];
        double dDT = dD[coordIdx];
        double dET = dE[coordIdx];
        double dFT = dF[coordIdx];

        double dA1 = 2 * (A*dET) - (dDT*F + D*dFT);
        double dA2 = 2 * (B*dDT) - 2 * (A*dET) + ((dDT - dET)*F + (D - E)*dFT);
        double dA3 = - 2 * (B*dDT) + (dFT*(E + F) + F*(dET + dFT));

        double dB1 = 4 * (A*dCT) - 2 * D*dDT;
        double dB2 = 4 * (A*dCT) - 2 * D*dDT + 4 * (B*dCT) - 2 * E*dET + 4 * (dCT*F + C*dFT) - 2 * (dDT*E + D*dET);
        double dB3 = - 2 * F*dFT + 4 * (B*dCT) - 2 * E*dET + 4 * (B*dDT) - 2 * (dET*F + E*dFT);
        double dBB1 = dB1 / 2 / BB1;
        double dBB2 = dB2 / 2 / BB2;
        double dBB3 = dB3 / 2 / BB3;

        double dC1 = dDT;
        double dC2 = dDT + dET + 2 * dFT;
        double dC3 = dET + dFT;
        double dD1 = dDT;
        double dD2 = dDT + dET;
        double dD3 = dET + dFT;

        double dE1 = (BB1*dC1 - C1*dBB1) / (B1 + C1*C1);
        double dE2 = (BB2*dC2 - C2*dBB2) / (B2 + C2*C2);
        double dE3 = (BB3*dC3 - C3*dBB3) / (B3 + C3*C3);
        double dF1 = (BB1*dD1 - D1*dBB1) / (B1 + D1*D1);
        double dF2 = (BB2*dD2 - D2*dBB2) / (B2 + D2*D2);
        double dF3 = (BB3*dD3 - D3*dBB3) / (B3 + D3*D3);

        double dG1 = (BB1*dA1 - A1*dBB1) / B1;
        double dG2 = (BB2*dA2 - A2*dBB2) / B2;
        double dG3 = (BB3*dA3 - A3*dBB3) / B3;

        double dNumerator = dG1*(E1 - F1) + G1*(dE1 - dF1) + dG2*(E2 - F2) + G2*(dE2 - dF2) + dG3*(E3 - F3) + G3*(dE3 - dF3);
        double dDenominator = B * 2 * D*dDT + A*(-4 * (B*dCT) + 2 * E*dET) + dFT*(-D*E + C*F) + F*(-(dDT*E + D*dET) + dCT*F + C*dFT);

        dH[coordIdx] = (denominator*dNumerator - numerator*dDenominator) / (denominator*denominator);
    }

    Vec3 forceBead = (-kExVol * 2 * area) * dH;

    return forceBead;

}
