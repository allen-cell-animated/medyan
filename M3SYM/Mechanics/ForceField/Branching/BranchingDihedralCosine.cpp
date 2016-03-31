
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "BranchingDihedralCosine.h"

#include "Bead.h"

#include "MathFunctions.h"

using namespace mathfunc;

double BranchingDihedralCosine::energy(Bead* b1, Bead* b2, Bead* b3, Bead* b4,
                                       double kDihed, double position){
    
    
    vector<double> n1 = vectorProduct(midPointCoordinate(b1->coordinate, b2->coordinate, position), b2->coordinate,
                                      midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate);
    
    
    vector<double> n2 = vectorProduct(b3->coordinate, b4->coordinate,
                        midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate);
    
    auto n1_norm = normalizedVector(n1);
    auto n2_norm = normalizedVector(n2);
    
    double n1n2 = dotProduct(n1_norm, n2_norm);
    
    return kDihed * ( 1 - n1n2 );
}


double BranchingDihedralCosine::energy(Bead* b1, Bead* b2, Bead* b3, Bead* b4,
                                       double kDihed, double position, double d){
    
    vector<double> zero (3,0); //Aux zero vector;
    
    vector<double> n1 = vectorProductStretched(
    midPointCoordinateStretched(b1->coordinate, b1->force, b2->coordinate, b2->force, position, d), zero,  b2->coordinate, b2->force,
    midPointCoordinateStretched(b1->coordinate, b1->force, b2->coordinate, b2->force, position, d), zero, b3->coordinate, b3->force, d);
    
    vector<double> n2 = vectorProductStretched(b3->coordinate,b3->force, b4->coordinate, b4->force,
    midPointCoordinateStretched(b1->coordinate, b1->force, b2->coordinate, b2->force, position, d), zero, b3->coordinate, b3->force, d);
    
    auto n1_norm = normalizedVector(n1);
    auto n2_norm = normalizedVector(n2);
    
    double n1n2 = dotProduct(n1_norm, n2_norm);
    
    return kDihed * ( 1 - n1n2 );
}

void BranchingDihedralCosine::forces(Bead* b1, Bead* b2, Bead* b3, Bead* b4,
                                     double kDihed, double position){
    
    vector<double> n1 = vectorProduct(midPointCoordinate(b1->coordinate, b2->coordinate, position), b2->coordinate,
                                      midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate);
    
    
    vector<double> n2 = vectorProduct(b3->coordinate, b4->coordinate,
                                      midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate);
    
    double N1 = sqrt(dotProduct(n1, n1));
    double N2 = sqrt(dotProduct(n2, n2));
    double n1n2 = dotProduct(n1, n2);
    
    double f0 = kDihed/N1/N2;
    
    double NN1 = n1n2/N1/N1;
    double NN2 = n1n2/N2/N2;
    
    double X = sqrt(scalarProduct(midPointCoordinate(b1->coordinate, b2->coordinate, position), b2->coordinate,
                                  midPointCoordinate(b1->coordinate, b2->coordinate, position), b2->coordinate));
    
    double D = sqrt(scalarProduct(midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate,
                                  midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate));
    
    double Y = sqrt(scalarProduct(b3->coordinate, b4->coordinate,
                                  b3->coordinate, b4->coordinate));
    
    vector<double> zero (3,0); //Aux zero vector;
    
    double n2x = scalarProduct(zero, n2, midPointCoordinate(b1->coordinate, b2->coordinate, position), b2->coordinate);
    
    double n1y = scalarProduct(zero, n1, b3->coordinate, b4->coordinate);
    
    double xd = scalarProduct(midPointCoordinate(b1->coordinate, b2->coordinate, position), b2->coordinate,
                              midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate);
    
    double yd = scalarProduct(b3->coordinate, b4->coordinate,
                              midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate);
    
    double xx = scalarProduct(midPointCoordinate(b1->coordinate, b2->coordinate, position), b2->coordinate,
                              midPointCoordinate(b1->coordinate, b2->coordinate, position), b2->coordinate);
    
    double xy = scalarProduct(midPointCoordinate(b1->coordinate, b2->coordinate, position),
                              b2->coordinate, b3->coordinate, b4->coordinate);
    
    double yy = scalarProduct(b3->coordinate, b4->coordinate, b3->coordinate, b4->coordinate);
    
    double XD = n2x/D/X/X/X;
    double X1 = -NN2*xd/D/X + yd/D/Y + yd/D/D/X/Y;
    double X2 = xd*yd/D/D/X/X/X/Y;
    double Y1 = -xd/D/X - xd/D/D/X/Y + NN1*yd/D/Y;
    double Y2 = xd*yd/D/D/X/Y/Y/Y;
    double D1 = NN2*xx/D/X - xy/D/X-xy/D/Y - 2*xy/D/D/X/Y + NN1*yy/D/Y;
    double D2 = xd*xy/D/D/X/X/X/Y;
    double YD = n1y/D/Y/Y/Y;
    
    //force on b1:
    b1->force[0] += f0*(- (1 - position)*XD*(1-position)*( (b2->coordinate[1] - b1->coordinate[1])*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) - (b2->coordinate[2] - b1->coordinate[2])*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) ) + (1 - position)*(X1 - X2)*(1-position)*(b2->coordinate[0] - b1->coordinate[0]) - (1 - position)*Y1*(b4->coordinate[0] - b3->coordinate[0]) + (1 - position)*(D1 + D2)*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]));
    
    
    b1->force[1] += f0*(- (1 - position)*XD*(1-position)*( (b2->coordinate[2] - b1->coordinate[2])*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) - (b2->coordinate[0] - b1->coordinate[0])*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) ) + (1 - position)*(X1 - X2)*(1-position)*(b2->coordinate[1] - b1->coordinate[1]) - (1 - position)*Y1*(b4->coordinate[1] - b3->coordinate[1]) + (1 - position)*(D1 + D2)*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]));
    
    b1->force[2] += f0*(- (1 - position)*XD*(1-position)*( (b2->coordinate[0] - b1->coordinate[0])*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) - (b2->coordinate[1] - b1->coordinate[1])*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) ) + (1 - position)*(X1 - X2)*(1-position)*(b2->coordinate[2] - b1->coordinate[2]) - (1 - position)*Y1*(b4->coordinate[2] - b3->coordinate[2]) + (1 - position)*(D1 + D2)*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]));
    
    
    //force on b2:
    b2->force[0] += f0*( (1 - position)*XD*(1-position)*( (b2->coordinate[1] - b1->coordinate[1])*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) - (b2->coordinate[2] - b1->coordinate[2])*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) ) + (X2 + position*(X1 - X2))*(1-position)*(b2->coordinate[0] - b1->coordinate[0]) - position*Y1*(b4->coordinate[0] - b3->coordinate[0]) + (position*(D1 + D2) - D2)*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) );
    
    b2->force[1] += f0*( (1 - position)*XD*(1-position)*( (b2->coordinate[2] - b1->coordinate[2])*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) - (b2->coordinate[0] - b1->coordinate[0])*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) ) + (X2 + position*(X1 - X2))*(1-position)*(b2->coordinate[1] - b1->coordinate[1]) - position*Y1*(b4->coordinate[1] - b3->coordinate[1]) + (position*(D1 + D2) - D2)*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) );
    
    b2->force[2] += f0*( (1 - position)*XD*(1-position)*( (b2->coordinate[0] - b1->coordinate[0])*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) - (b2->coordinate[1] - b1->coordinate[1])*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) ) + (X2 + position*(X1 - X2))*(1-position)*(b2->coordinate[2] - b1->coordinate[2]) - position*Y1*(b4->coordinate[2] - b3->coordinate[2]) + (position*(D1 + D2) - D2)*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) );
    
    //force on b3:
    b3->force[0] += f0*(-YD*( (b4->coordinate[1] - b3->coordinate[1])*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) - (b4->coordinate[2] - b3->coordinate[2])*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) ) - X1*(1-position)*(b2->coordinate[0] - b1->coordinate[0]) + (Y1 - Y2)*(b4->coordinate[0] - b3->coordinate[0]) + (D2 - D1)*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]));
    
    b3->force[1] += f0*(-YD*( (b4->coordinate[2] - b3->coordinate[2])*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) - (b4->coordinate[0] - b3->coordinate[0])*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) ) - X1*(1-position)*(b2->coordinate[1] - b1->coordinate[1]) + (Y1 - Y2)*(b4->coordinate[1] - b3->coordinate[1]) + (D2 - D1)*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]));
    
    b3->force[2] += f0*(-YD*( (b4->coordinate[0] - b3->coordinate[0])*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) - (b4->coordinate[1] - b3->coordinate[1])*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) ) - X1*(1-position)*(b2->coordinate[2] - b1->coordinate[2]) + (Y1 - Y2)*(b4->coordinate[2] - b3->coordinate[2]) + (D2 - D1)*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]));
    
    
    //force on b4:
    b4->force[0] +=f0*( YD*( (b4->coordinate[1] - b3->coordinate[1])*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) - (b4->coordinate[2] - b3->coordinate[2])*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) ) + Y2*(b4->coordinate[0] - b3->coordinate[0]) - D2*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) );
    
    b4->force[1] +=f0*( YD*( (b4->coordinate[2] - b3->coordinate[2])*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) - (b4->coordinate[0] - b3->coordinate[0])*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) ) + Y2*(b4->coordinate[1] - b3->coordinate[1]) - D2*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) );
    
    b4->force[2] +=f0*( YD*( (b4->coordinate[0] - b3->coordinate[0])*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) - (b4->coordinate[1] - b3->coordinate[1])*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) ) + Y2*(b4->coordinate[2] - b3->coordinate[2]) - D2*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) );

    
}

void BranchingDihedralCosine::forcesAux(Bead* b1, Bead* b2, Bead* b3, Bead* b4,
                                     double kDihed, double position){
    
    
    vector<double> n1 = vectorProduct(midPointCoordinate(b1->coordinate, b2->coordinate, position), b2->coordinate,
                                      midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate);
    
    
    vector<double> n2 = vectorProduct(b3->coordinate, b4->coordinate,
                                      midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate);
    
    double N1 = sqrt(dotProduct(n1, n1));
    double N2 = sqrt(dotProduct(n2, n2));
    double n1n2 = dotProduct(n1, n2);
    
    double f0 = kDihed/N1/N2;
    
    double NN1 = n1n2/N1/N1;
    double NN2 = n1n2/N2/N2;
    
    double X = sqrt(scalarProduct(midPointCoordinate(b1->coordinate, b2->coordinate, position), b2->coordinate,
                                  midPointCoordinate(b1->coordinate, b2->coordinate, position), b2->coordinate));
    
    double D = sqrt(scalarProduct(midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate,
                                  midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate));
    
    double Y = sqrt(scalarProduct(b3->coordinate, b4->coordinate,
                                  b3->coordinate, b4->coordinate));
    
    vector<double> zero (3,0); //Aux zero vector;
    
    double n2x = scalarProduct(zero, n2, midPointCoordinate(b1->coordinate, b2->coordinate, position), b2->coordinate);
    
    double n1y = scalarProduct(zero, n1, b3->coordinate, b4->coordinate);
    
    double xd = scalarProduct(midPointCoordinate(b1->coordinate, b2->coordinate, position), b2->coordinate,
                              midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate);
    
    double yd = scalarProduct(b3->coordinate, b4->coordinate,
                              midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate);
    
    double xx = scalarProduct(midPointCoordinate(b1->coordinate, b2->coordinate, position),
                              b2->coordinate, midPointCoordinate(b1->coordinate, b2->coordinate, position),
                              b2->coordinate);
    
    double xy = scalarProduct(midPointCoordinate(b1->coordinate, b2->coordinate, position),
                              b2->coordinate, b3->coordinate, b4->coordinate);
    
    double yy = scalarProduct(b3->coordinate, b4->coordinate, b3->coordinate, b4->coordinate);
    
    double XD = n2x/D/X/X/X;
    double X1 = -NN2*xd/D/X + yd/D/Y + yd/D/D/X/Y;
    double X2 = xd*yd/D/D/X/X/X/Y;
    double Y1 = -xd/D/X - xd/D/D/X/Y + NN1*yd/D/Y;
    double Y2 = xd*yd/D/D/X/Y/Y/Y;
    double D1 = NN2*xx/D/X - xy/D/X-xy/D/Y - 2*xy/D/D/X/Y + NN1*yy/D/Y;
    double D2 = xd*xy/D/D/X/X/X/Y;
    double YD = n1y/D/Y/Y/Y;
    
    //forceAux on b1:
    b1->forceAux[0] += f0*(- (1 - position)*XD*(1-position)*( (b2->coordinate[1] - b1->coordinate[1])*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) - (b2->coordinate[2] - b1->coordinate[2])*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) ) + (1 - position)*(X1 - X2)*(1-position)*(b2->coordinate[0] - b1->coordinate[0]) - (1 - position)*Y1*(b4->coordinate[0] - b3->coordinate[0]) + (1 - position)*(D1 + D2)*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]));
    
    b1->forceAux[1] += f0*(- (1 - position)*XD*(1-position)*( (b2->coordinate[2] - b1->coordinate[2])*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) - (b2->coordinate[0] - b1->coordinate[0])*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) ) + (1 - position)*(X1 - X2)*(1-position)*(b2->coordinate[1] - b1->coordinate[1]) - (1 - position)*Y1*(b4->coordinate[1] - b3->coordinate[1]) + (1 - position)*(D1 + D2)*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]));
    
    b1->forceAux[2] += f0*(- (1 - position)*XD*(1-position)*( (b2->coordinate[0] - b1->coordinate[0])*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) - (b2->coordinate[1] - b1->coordinate[1])*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) ) + (1 - position)*(X1 - X2)*(1-position)*(b2->coordinate[2] - b1->coordinate[2]) - (1 - position)*Y1*(b4->coordinate[2] - b3->coordinate[2]) + (1 - position)*(D1 + D2)*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]));
    
    
    //forceAux on b2:
    b2->forceAux[0] += f0*( (1 - position)*XD*(1-position)*( (b2->coordinate[1] - b1->coordinate[1])*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) - (b2->coordinate[2] - b1->coordinate[2])*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) ) + (X2 + position*(X1 - X2))*(1-position)*(b2->coordinate[0] - b1->coordinate[0]) - position*Y1*(b4->coordinate[0] - b3->coordinate[0]) + (position*(D1 + D2) - D2)*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) );
    
    b2->forceAux[1] += f0*( (1 - position)*XD*(1-position)*( (b2->coordinate[2] - b1->coordinate[2])*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) - (b2->coordinate[0] - b1->coordinate[0])*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) ) + (X2 + position*(X1 - X2))*(1-position)*(b2->coordinate[1] - b1->coordinate[1]) - position*Y1*(b4->coordinate[1] - b3->coordinate[1]) + (position*(D1 + D2) - D2)*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) );
    
    b2->forceAux[2] += f0*( (1 - position)*XD*(1-position)*( (b2->coordinate[0] - b1->coordinate[0])*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) - (b2->coordinate[1] - b1->coordinate[1])*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) ) + (X2 + position*(X1 - X2))*(1-position)*(b2->coordinate[2] - b1->coordinate[2]) - position*Y1*(b4->coordinate[2] - b3->coordinate[2]) + (position*(D1 + D2) - D2)*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) );
    
    //forceAux on b3:
    b3->forceAux[0] += f0*(-YD*( (b4->coordinate[1] - b3->coordinate[1])*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) - (b4->coordinate[2] - b3->coordinate[2])*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) ) - X1*(1-position)*(b2->coordinate[0] - b1->coordinate[0]) + (Y1 - Y2)*(b4->coordinate[0] - b3->coordinate[0]) + (D2 - D1)*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]));
    
    b3->forceAux[1] += f0*(-YD*( (b4->coordinate[2] - b3->coordinate[2])*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) - (b4->coordinate[0] - b3->coordinate[0])*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) ) - X1*(1-position)*(b2->coordinate[1] - b1->coordinate[1]) + (Y1 - Y2)*(b4->coordinate[1] - b3->coordinate[1]) + (D2 - D1)*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]));
    
    b3->forceAux[2] += f0*(-YD*( (b4->coordinate[0] - b3->coordinate[0])*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) - (b4->coordinate[1] - b3->coordinate[1])*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) ) - X1*(1-position)*(b2->coordinate[2] - b1->coordinate[2]) + (Y1 - Y2)*(b4->coordinate[2] - b3->coordinate[2]) + (D2 - D1)*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]));
    
    
    //forceAux on b4:
    b4->forceAux[0] +=f0*( YD*( (b4->coordinate[1] - b3->coordinate[1])*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) - (b4->coordinate[2] - b3->coordinate[2])*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) ) + Y2*(b4->coordinate[0] - b3->coordinate[0]) - D2*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) );
    
    b4->forceAux[1] +=f0*( YD*( (b4->coordinate[2] - b3->coordinate[2])*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) - (b4->coordinate[0] - b3->coordinate[0])*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) ) + Y2*(b4->coordinate[1] - b3->coordinate[1]) - D2*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) );
    
    b4->forceAux[2] +=f0*( YD*( (b4->coordinate[0] - b3->coordinate[0])*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) - (b4->coordinate[1] - b3->coordinate[1])*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) ) + Y2*(b4->coordinate[2] - b3->coordinate[2]) - D2*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) );
}