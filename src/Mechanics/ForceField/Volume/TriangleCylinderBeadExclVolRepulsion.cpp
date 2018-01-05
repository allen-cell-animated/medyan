
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

// Note: These functions require that the area of the triangle has already been calculated

double TriangleCylinderBeadExclVolRepulsion::energy(Triangle* t, Bead* b,
                                                    double kExVol) {
    
    auto& v = t->getVertices();

    auto c0 = v[0]->coordinate;
    auto c1 = v[1]->coordinate;
    auto c2 = v[2]->coordinate;
    auto cb = b->coordinate;
        
    //check if in same plane
    if(areInPlane(c0, c1, c2, cb)) {
        
        // slightly move point. Using negative number here to move "into the cell".
        cb = movePointOutOfPlane(c0, c1, c2, cb, 4, -0.01);
    }
    
    double A = scalarProduct(c0, c1, c0, c1);
    double B = scalarProduct(c1, c2, c1, c2);
    double C = scalarProduct(cb, c0, cb, c0);
    double D = 2 * scalarProduct(c0, c1, cb, c0);
    double E = 2 * scalarProduct(c1, c2, cb, c0);
    double F = 2 * scalarProduct(c0, c1, c1, c2);

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
    
    double energy = kExVol * t->getGTriangle()->getArea() * H;
    
    return energy;
}

double TriangleCylinderBeadExclVolRepulsion::energy(Triangle* t, Bead* b,
                                                    double kExVol, double d) {
    
    auto& v = t->getVertices();

    auto& c0Raw = v[0]->coordinate;
    auto& c1Raw = v[1]->coordinate;
    auto& c2Raw = v[2]->coordinate;
    auto& cbRaw = b->coordinate;

    auto& f0 = v[0]->force;
    auto& f1 = v[1]->force;
    auto& f2 = v[2]->force;
    auto& fb = b->force;

    vector<double> c0 = {c0Raw[0] + d * f0[0], c0Raw[1] + d * f0[1], c0Raw[2] + d * f0[2]};
    vector<double> c1 = {c1Raw[0] + d * f1[0], c1Raw[1] + d * f1[1], c1Raw[2] + d * f1[2]};
    vector<double> c2 = {c2Raw[0] + d * f2[0], c2Raw[1] + d * f2[1], c2Raw[2] + d * f2[2]};
    vector<double> cb = {cbRaw[0] + d * fb[0], cbRaw[1] + d * fb[1], cbRaw[2] + d * fb[2]};

    //check if in same plane
    if(areInPlane(c0, c1, c2, cb)) {
        
        // slightly move point. Using negative number here to move "into the cell".
        cb = movePointOutOfPlane(c0, c1, c2, cb, 4, -0.01);
    }
    
    double A = scalarProduct(c0, c1, c0, c1);
    double B = scalarProduct(c1, c2, c1, c2);
    double C = scalarProduct(cb, c0, cb, c0);
    double D = 2 * scalarProduct(c0, c1, cb, c0);
    double E = 2 * scalarProduct(c1, c2, cb, c0);
    double F = 2 * scalarProduct(c0, c1, c1, c2);

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
    
    double energy = kExVol * t->getGTriangle()->getStretchedArea() * H;
    
    return energy;
}

void TriangleCylinderBeadExclVolRepulsion::forces(Triangle* t, Bead* b,
                                                  double kExVol) {
    
    auto& v = t->getVertices();

    auto c0 = v[0]->coordinate;
    auto c1 = v[1]->coordinate;
    auto c2 = v[2]->coordinate;
    auto cb = b->coordinate;
        
    //check if in same plane
    if(areInPlane(c0, c1, c2, cb)) {
        
        // slightly move point. Using negative number here to move "into the cell".
        cb = movePointOutOfPlane(c0, c1, c2, cb, 4, -0.01);
    }
    
    double A = scalarProduct(c0, c1, c0, c1);
    double B = scalarProduct(c1, c2, c1, c2);
    double C = scalarProduct(cb, c0, cb, c0);
    double D = 2 * scalarProduct(c0, c1, cb, c0);
    double E = 2 * scalarProduct(c1, c2, cb, c0);
    double F = 2 * scalarProduct(c0, c1, c1, c2);

    // Derivative index is 0, 1, 2, b
    array<array<double, 3>, 4> dA = {{
        {{ 2 * (c0[0] - c1[0]), 2 * (c0[1] - c1[1]), 2 * (c0[2] - c1[2]) }},
        {{ 2 * (c1[0] - c0[0]), 2 * (c1[1] - c0[1]), 2 * (c1[2] - c0[2]) }}
    }};
    array<array<double, 3>, 4> dB = {{
        {{ 0, 0, 0 }},
        {{ 2 * (c1[0] - c2[0]), 2 * (c1[1] - c2[1]), 2 * (c1[2] - c2[2]) }},
        {{ 2 * (c2[0] - c1[0]), 2 * (c2[1] - c1[1]), 2 * (c2[2] - c1[2]) }}
    }};
    array<array<double, 3>, 4> dC = {{
        {{ 2 * (c0[0] - cb[0]), 2 * (c0[1] - cb[1]), 2 * (c0[2] - cb[2]) }},
        {{ 0, 0, 0 }},
        {{ 0, 0, 0 }},
        {{ 2 * (cb[0] - c0[0]), 2 * (cb[1] - c0[1]), 2 * (cb[2] - c0[2]) }}
    }};
    array<array<double, 3>, 4> dD = {{
        {{ 2 * (c1[0] + cb[0] - 2*c0[0]), 2 * (c1[1] + cb[1] - 2*c0[1]), 2 * (c1[2] + cb[2] - 2*c0[2]) }},
        {{ 2 * (c0[0] - cb[0]), 2 * (c0[1] - cb[1]), 2 * (c0[2] - cb[2]) }},
        {{ 0, 0, 0 }},
        {{ 2 * (c0[0] - c1[0]), 2 * (c0[1] - c1[1]), 2 * (c0[2] - c1[2]) }}
    }};
    array<array<double, 3>, 4> dE = {{
        {{ 2 * (c2[0] - c1[0]), 2 * (c2[1] - c1[1]), 2 * (c2[2] - c1[2]) }},
        {{ 2 * (cb[0] - c0[0]), 2 * (cb[1] - c0[1]), 2 * (cb[2] - c0[2]) }},
        {{ 2 * (c0[0] - cb[0]), 2 * (c0[1] - cb[1]), 2 * (c0[2] - cb[2]) }},
        {{ 2 * (c1[0] - c2[0]), 2 * (c1[1] - c2[1]), 2 * (c1[2] - c2[2]) }}
    }};
    array<array<double, 3>, 4> dF = {{
        {{ 2 * (c1[0] - c2[0]), 2 * (c1[1] - c2[1]), 2 * (c1[2] - c2[2]) }},
        {{ 2 * (c0[0] + c2[0] - 2*c1[0]), 2 * (c0[1] + c2[1] - 2*c1[1]), 2 * (c0[2] + c2[2] - 2*c1[2]) }},
        {{ 2 * (c1[0] - c0[0]), 2 * (c1[1] - c0[1]), 2 * (c1[2] - c0[2]) }}
    }};

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
    
    array<array<double, 3>, 4> dH = {};
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

    double area = t->getGTriangle()->getArea();
    auto& dArea = t->getGTriangle()->getDArea();
    
    for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
        b->force[coordIdx] -= kExVol * dH[3][coordIdx] * area;
        for(size_t vtxIdx = 0; vtxIdx < 3; ++vtxIdx) {
            v[vtxIdx]->force[coordIdx] -= kExVol * (dArea[vtxIdx][coordIdx]*H + dH[vtxIdx][coordIdx]*area);
        }
    }
}

void TriangleCylinderBeadExclVolRepulsion::forcesAux(Triangle* t, Bead* b,
                                                     double kExVol) {
    
    auto& v = t->getVertices();

    auto c0 = v[0]->coordinate;
    auto c1 = v[1]->coordinate;
    auto c2 = v[2]->coordinate;
    auto cb = b->coordinate;
        
    //check if in same plane
    if(areInPlane(c0, c1, c2, cb)) {
        
        // slightly move point. Using negative number here to move "into the cell".
        cb = movePointOutOfPlane(c0, c1, c2, cb, 4, -0.01);
    }
    
    double A = scalarProduct(c0, c1, c0, c1);
    double B = scalarProduct(c1, c2, c1, c2);
    double C = scalarProduct(cb, c0, cb, c0);
    double D = 2 * scalarProduct(c0, c1, cb, c0);
    double E = 2 * scalarProduct(c1, c2, cb, c0);
    double F = 2 * scalarProduct(c0, c1, c1, c2);

    // Derivative index is 0, 1, 2, b
    array<array<double, 3>, 4> dA = {{
        {{ 2 * (c0[0] - c1[0]), 2 * (c0[1] - c1[1]), 2 * (c0[2] - c1[2]) }},
        {{ 2 * (c1[0] - c0[0]), 2 * (c1[1] - c0[1]), 2 * (c1[2] - c0[2]) }}
    }};
    array<array<double, 3>, 4> dB = {{
        {{ 0, 0, 0 }},
        {{ 2 * (c1[0] - c2[0]), 2 * (c1[1] - c2[1]), 2 * (c1[2] - c2[2]) }},
        {{ 2 * (c2[0] - c1[0]), 2 * (c2[1] - c1[1]), 2 * (c2[2] - c1[2]) }}
    }};
    array<array<double, 3>, 4> dC = {{
        {{ 2 * (c0[0] - cb[0]), 2 * (c0[1] - cb[1]), 2 * (c0[2] - cb[2]) }},
        {{ 0, 0, 0 }},
        {{ 0, 0, 0 }},
        {{ 2 * (cb[0] - c0[0]), 2 * (cb[1] - c0[1]), 2 * (cb[2] - c0[2]) }}
    }};
    array<array<double, 3>, 4> dD = {{
        {{ 2 * (c1[0] + cb[0] - 2*c0[0]), 2 * (c1[1] + cb[1] - 2*c0[1]), 2 * (c1[2] + cb[2] - 2*c0[2]) }},
        {{ 2 * (c0[0] - cb[0]), 2 * (c0[1] - cb[1]), 2 * (c0[2] - cb[2]) }},
        {{ 0, 0, 0 }},
        {{ 2 * (c0[0] - c1[0]), 2 * (c0[1] - c1[1]), 2 * (c0[2] - c1[2]) }}
    }};
    array<array<double, 3>, 4> dE = {{
        {{ 2 * (c2[0] - c1[0]), 2 * (c2[1] - c1[1]), 2 * (c2[2] - c1[2]) }},
        {{ 2 * (cb[0] - c0[0]), 2 * (cb[1] - c0[1]), 2 * (cb[2] - c0[2]) }},
        {{ 2 * (c0[0] - cb[0]), 2 * (c0[1] - cb[1]), 2 * (c0[2] - cb[2]) }},
        {{ 2 * (c1[0] - c2[0]), 2 * (c1[1] - c2[1]), 2 * (c1[2] - c2[2]) }}
    }};
    array<array<double, 3>, 4> dF = {{
        {{ 2 * (c1[0] - c2[0]), 2 * (c1[1] - c2[1]), 2 * (c1[2] - c2[2]) }},
        {{ 2 * (c0[0] + c2[0] - 2*c1[0]), 2 * (c0[1] + c2[1] - 2*c1[1]), 2 * (c0[2] + c2[2] - 2*c1[2]) }},
        {{ 2 * (c1[0] - c0[0]), 2 * (c1[1] - c0[1]), 2 * (c1[2] - c0[2]) }}
    }};

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
    
    array<array<double, 3>, 4> dH = {};
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

    double area = t->getGTriangle()->getArea();
    auto& dArea = t->getGTriangle()->getDArea();
    
    for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
        b->forceAux[coordIdx] -= kExVol * dH[3][coordIdx] * area;
        for(size_t vtxIdx = 0; vtxIdx < 3; ++vtxIdx) {
            v[vtxIdx]->forceAux[coordIdx] -= kExVol * (dArea[vtxIdx][coordIdx]*H + dH[vtxIdx][coordIdx]*area);
        }
    }
}

