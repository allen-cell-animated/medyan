
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "BranchingDihedralCosine.h"

#include "Bead.h"

#include "MathFunctions.h"

using namespace mathfunc;

double BranchingDihedralCosine::energy(Bead* b1, Bead* b2, Bead* b3, Bead* b4,
                                       double kDihed, double position){
    
    
    vector<double> n1 = vectorProduct(midPointCoordinate(b1->coordinate, b2->coordinate, position),b2->coordinate,
                                      midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate);
    
    vector<double> n2 = vectorProduct(b3->coordinate, b4->coordinate,
                                      midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate);
    
    double norm_x = sqrt(scalarProduct(midPointCoordinate(b1->coordinate, b2->coordinate, position),b2->coordinate,
                                   midPointCoordinate(b1->coordinate, b2->coordinate, position),b2->coordinate));
    double norm_d = sqrt(scalarProduct(midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate,
                                  midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate));
    
    double norm_y = sqrt(scalarProduct(b3->coordinate, b4->coordinate,
                                      b3->coordinate, b4->coordinate));
    
    double N1 = 1/norm_x/norm_d;
    double N2 = 1/norm_y/norm_d;
    
    double n1n2 = dotProduct(n1, n2);
   
    
    return kDihed * ( 1 - n1n2*N1*N2 );
    
}


double BranchingDihedralCosine::energy(Bead* b1, Bead* b2, Bead* b3, Bead* b4,
                                       double kDihed, double position, double d){
    
    vector<double> zero (3,0); //Aux zero vector;
    
    vector<double> n1 = vectorProductStretched(midPointCoordinateStretched(b1->coordinate,b1->force, b2->coordinate,b2->force, position, d),zero ,b2->coordinate, b2->force,
                                      midPointCoordinateStretched(b1->coordinate,b1->force, b2->coordinate,b2->force, position, d),zero, b3->coordinate, b3->force, d);
    
    vector<double> n2 = vectorProductStretched(b3->coordinate,b3->force, b4->coordinate,b4->force,
                                      midPointCoordinateStretched(b1->coordinate,b1->force, b2->coordinate,b2->force, position, d),zero, b3->coordinate,b3->force, d);
    
    double norm_x = sqrt(scalarProductStretched(midPointCoordinateStretched(b1->coordinate,b1->force, b2->coordinate,b2->force, position, d),zero ,b2->coordinate, b2->force,
                                       midPointCoordinateStretched(b1->coordinate,b1->force, b2->coordinate,b2->force, position, d),zero ,b2->coordinate, b2->force,  d));
    
    double norm_d = sqrt(scalarProductStretched(midPointCoordinateStretched(b1->coordinate,b1->force, b2->coordinate,b2->force, position, d),zero ,b3->coordinate, b3->force,
                                                midPointCoordinateStretched(b1->coordinate,b1->force, b2->coordinate,b2->force, position, d),zero ,b3->coordinate, b3->force,d));
    
    double norm_y = sqrt(scalarProductStretched(b3->coordinate,b3->force, b4->coordinate,b4->force,
                                       b3->coordinate,b3->force, b4->coordinate,b4->force, d));
    
    double N1 = 1/norm_x/norm_d;
    double N2 = 1/norm_y/norm_d;
    
    double n1n2 = dotProduct(n1, n2);
    
    
    return kDihed * ( 1 - n1n2*N1*N2 );
}
    
void BranchingDihedralCosine::forces(Bead* b1, Bead* b2, Bead* b3, Bead* b4,
                                            double kDihed, double position){

    double a = scalarProduct(midPointCoordinate(b1->coordinate, b2->coordinate, position),b2->coordinate,
                             b3->coordinate, b4->coordinate);
    
    double b = scalarProduct(midPointCoordinate(b1->coordinate, b2->coordinate, position),b2->coordinate,
                             midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate);

    double c = scalarProduct(b3->coordinate, b4->coordinate,
                             midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate);
    
    
    vector<double> n1 = vectorProduct(midPointCoordinate(b1->coordinate, b2->coordinate, position),b2->coordinate,
                                      midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate);
    
    vector<double> n2 = vectorProduct(b3->coordinate, b4->coordinate,
                                      midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate);
    
    double n = dotProduct(n1, n2);
    
    
    double norm_x = sqrt(scalarProduct(midPointCoordinate(b1->coordinate, b2->coordinate, position),b2->coordinate,
                                       midPointCoordinate(b1->coordinate, b2->coordinate, position),b2->coordinate));
    double norm_d = sqrt(scalarProduct(midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate,
                                       midPointCoordinate(b1->coordinate, b2->coordinate, position), b3->coordinate));
    
    double norm_y = sqrt(scalarProduct(b3->coordinate, b4->coordinate,
                                       b3->coordinate, b4->coordinate));
    
    double N1 = 1/norm_x/norm_d;
    double N2 = 1/norm_y/norm_d;
    
    
    double X = norm_x/norm_d*n/N1;
    double Y = norm_y/norm_d*n/N2;
    
    double DX = norm_d/norm_x*n/N1;
    double DY = norm_d/norm_y*n/N2;
 
//
    
    b1->force[0] += kDihed*( (1-position)*N2*N2*((Y-a)*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) + b*(b4->coordinate[0] - b3->coordinate[0])) + (1-position)*N1*N1*( (X-a)*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) + (c+DX)*(1-position)*(b2->coordinate[0] - b1->coordinate[0]) ) );
    
    b1->force[1] += kDihed*( (1-position)*N2*N2*((Y-a)*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) + b*(b4->coordinate[1] - b3->coordinate[1])) + (1-position)*N1*N1*( (X-a)*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) + (c+DX)*(1-position)*(b2->coordinate[1] - b1->coordinate[1]) ) );
    
    b1->force[2] += kDihed*( (1-position)*N2*N2*((Y-a)*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) + b*(b4->coordinate[2] - b3->coordinate[2])) + (1-position)*N1*N1*( (X-a)*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) + (c+DX)*(1-position)*(b2->coordinate[2] - b1->coordinate[2]) ) );
    
 //
    
    b2->force[0] += kDihed*( position*N2*N2*((Y-a)*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) -b*(b4->coordinate[0] - b3->coordinate[0])) + position*N1*N1*((X-a) + c*(1-position)*(b2->coordinate[0] - b1->coordinate[0])) - (1-position)*N1*N1*DX*(1-position)*(b2->coordinate[0] - b1->coordinate[0]));
    
    b2->force[1] += kDihed*( position*N2*N2*((Y-a)*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) -b*(b4->coordinate[1] - b3->coordinate[1])) + position*N1*N1*((X-a) + c*(1-position)*(b2->coordinate[1] - b1->coordinate[1])) - (1-position)*N1*N1*DX*(1-position)*(b2->coordinate[1] - b1->coordinate[1]));
    
    b2->force[2] += kDihed*( position*N2*N2*((Y-a)*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) -b*(b4->coordinate[2] - b3->coordinate[2])) + position*N1*N1*((X-a) + c*(1-position)*(b2->coordinate[2] - b1->coordinate[2])) - (1-position)*N1*N1*DX*(1-position)*(b2->coordinate[2] - b1->coordinate[2]));
//
    
    
    b3->force[0] += kDihed*( N2*N2*( (DY-b)*(b4->coordinate[0] - b3->coordinate[0]) +(a-Y)*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) ) + N1*N1*( (a-X)*(b3->coordinate[0] - (1-position)*b1->coordinate[0] - position*b2->coordinate[0]) -c*(1-position)*(b2->coordinate[0] - b1->coordinate[0])) );
    
    b3->force[1] += kDihed*( N2*N2*( (DY-b)*(b4->coordinate[1] - b3->coordinate[1]) +(a-Y)*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) ) + N1*N1*( (a-X)*(b3->coordinate[1] - (1-position)*b1->coordinate[1] - position*b2->coordinate[1]) -c*(1-position)*(b2->coordinate[0] - b1->coordinate[1])) );
    
    b3->force[2] += kDihed*( N2*N2*( (DY-b)*(b4->coordinate[2] - b3->coordinate[2]) +(a-Y)*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) ) + N1*N1*( (a-X)*(b3->coordinate[2] - (1-position)*b1->coordinate[2] - position*b2->coordinate[2]) -c*(1-position)*(b2->coordinate[2] - b1->coordinate[2])) );
    
//
    
    b4->force[0] += -kDihed*DY*(b4->coordinate[0] - b3->coordinate[0]);
    
    b4->force[1] += -kDihed*DY*(b4->coordinate[1] - b3->coordinate[1]);
    
    b4->force[2] += -kDihed*DY*(b4->coordinate[2] - b3->coordinate[2]);
    
}

void BranchingDihedralCosine::forcesAux(Bead* b1, Bead* b2, Bead* b3, Bead* b4,
                                     double kDihed, double position){
    
    double a = scalarProduct(midPointCoordinate(b1->coordinateAux, b2->coordinateAux, position),b2->coordinateAux,
                             b3->coordinateAux, b4->coordinateAux);
    
    double b = scalarProduct(midPointCoordinate(b1->coordinateAux, b2->coordinateAux, position),b2->coordinateAux,
                             midPointCoordinate(b1->coordinateAux, b2->coordinateAux, position), b3->coordinateAux);
    
    double c = scalarProduct(b3->coordinateAux, b4->coordinateAux,
                             midPointCoordinate(b1->coordinateAux, b2->coordinateAux, position), b3->coordinateAux);
    
    
    vector<double> n1 = vectorProduct(midPointCoordinate(b1->coordinateAux, b2->coordinateAux, position),b2->coordinateAux,
                                      midPointCoordinate(b1->coordinateAux, b2->coordinateAux, position), b3->coordinateAux);
    
    vector<double> n2 = vectorProduct(b3->coordinateAux, b4->coordinateAux,
                                      midPointCoordinate(b1->coordinateAux, b2->coordinateAux, position), b3->coordinateAux);
    
    double n = dotProduct(n1, n2);
    
    
    double norm_x = sqrt(scalarProduct(midPointCoordinate(b1->coordinateAux, b2->coordinateAux, position),b2->coordinateAux,
                                       midPointCoordinate(b1->coordinateAux, b2->coordinateAux, position),b2->coordinateAux));
    double norm_d = sqrt(scalarProduct(midPointCoordinate(b1->coordinateAux, b2->coordinateAux, position), b3->coordinateAux,
                                       midPointCoordinate(b1->coordinateAux, b2->coordinateAux, position), b3->coordinateAux));
    
    double norm_y = sqrt(scalarProduct(b3->coordinateAux, b4->coordinateAux,
                                       b3->coordinateAux, b4->coordinateAux));
    
    double N1 = 1/norm_x/norm_d;
    double N2 = 1/norm_y/norm_d;
    
    
    double X = norm_x/norm_d*n/N1;
    double Y = norm_y/norm_d*n/N2;
    
    double DX = norm_d/norm_x*n/N1;
    double DY = norm_d/norm_y*n/N2;
    
    //
    
    b1->forceAux[0] += kDihed*( (1-position)*N2*N2*((Y-a)*(b3->coordinateAux[0] - (1-position)*b1->coordinateAux[0] - position*b2->coordinateAux[0]) + b*(b4->coordinateAux[0] - b3->coordinateAux[0])) + (1-position)*N1*N1*( (X-a)*(b3->coordinateAux[0] - (1-position)*b1->coordinateAux[0] - position*b2->coordinateAux[0]) + (c+DX)*(1-position)*(b2->coordinateAux[0] - b1->coordinateAux[0]) ) );
    
    b1->forceAux[1] += kDihed*( (1-position)*N2*N2*((Y-a)*(b3->coordinateAux[1] - (1-position)*b1->coordinateAux[1] - position*b2->coordinateAux[1]) + b*(b4->coordinateAux[1] - b3->coordinateAux[1])) + (1-position)*N1*N1*( (X-a)*(b3->coordinateAux[1] - (1-position)*b1->coordinateAux[1] - position*b2->coordinateAux[1]) + (c+DX)*(1-position)*(b2->coordinateAux[1] - b1->coordinateAux[1]) ) );
    
    b1->forceAux[2] += kDihed*( (1-position)*N2*N2*((Y-a)*(b3->coordinateAux[2] - (1-position)*b1->coordinateAux[2] - position*b2->coordinateAux[2]) + b*(b4->coordinateAux[2] - b3->coordinateAux[2])) + (1-position)*N1*N1*( (X-a)*(b3->coordinateAux[2] - (1-position)*b1->coordinateAux[2] - position*b2->coordinateAux[2]) + (c+DX)*(1-position)*(b2->coordinateAux[2] - b1->coordinateAux[2]) ) );
    
    //
    
    b2->forceAux[0] += kDihed*( position*N2*N2*((Y-a)*(b3->coordinateAux[0] - (1-position)*b1->coordinateAux[0] - position*b2->coordinateAux[0]) -b*(b4->coordinateAux[0] - b3->coordinateAux[0])) + position*N1*N1*((X-a) + c*(1-position)*(b2->coordinateAux[0] - b1->coordinateAux[0])) - (1-position)*N1*N1*DX*(1-position)*(b2->coordinateAux[0] - b1->coordinateAux[0]));
    
    b2->forceAux[1] += kDihed*( position*N2*N2*((Y-a)*(b3->coordinateAux[1] - (1-position)*b1->coordinateAux[1] - position*b2->coordinateAux[1]) -b*(b4->coordinateAux[1] - b3->coordinateAux[1])) + position*N1*N1*((X-a) + c*(1-position)*(b2->coordinateAux[1] - b1->coordinateAux[1])) - (1-position)*N1*N1*DX*(1-position)*(b2->coordinateAux[1] - b1->coordinateAux[1]));
    
    b2->forceAux[2] += kDihed*( position*N2*N2*((Y-a)*(b3->coordinateAux[2] - (1-position)*b1->coordinateAux[2] - position*b2->coordinateAux[2]) -b*(b4->coordinateAux[2] - b3->coordinateAux[2])) + position*N1*N1*((X-a) + c*(1-position)*(b2->coordinateAux[2] - b1->coordinateAux[2])) - (1-position)*N1*N1*DX*(1-position)*(b2->coordinateAux[2] - b1->coordinateAux[2]));
    //
    
    
    b3->forceAux[0] += kDihed*( N2*N2*( (DY-b)*(b4->coordinateAux[0] - b3->coordinateAux[0]) +(a-Y)*(b3->coordinateAux[0] - (1-position)*b1->coordinateAux[0] - position*b2->coordinateAux[0]) ) + N1*N1*( (a-X)*(b3->coordinateAux[0] - (1-position)*b1->coordinateAux[0] - position*b2->coordinateAux[0]) -c*(1-position)*(b2->coordinateAux[0] - b1->coordinateAux[0])) );
    
    b3->forceAux[1] += kDihed*( N2*N2*( (DY-b)*(b4->coordinateAux[1] - b3->coordinateAux[1]) +(a-Y)*(b3->coordinateAux[1] - (1-position)*b1->coordinateAux[1] - position*b2->coordinateAux[1]) ) + N1*N1*( (a-X)*(b3->coordinateAux[1] - (1-position)*b1->coordinateAux[1] - position*b2->coordinateAux[1]) -c*(1-position)*(b2->coordinateAux[0] - b1->coordinateAux[1])) );
    
    b3->forceAux[2] += kDihed*( N2*N2*( (DY-b)*(b4->coordinateAux[2] - b3->coordinateAux[2]) +(a-Y)*(b3->coordinateAux[2] - (1-position)*b1->coordinateAux[2] - position*b2->coordinateAux[2]) ) + N1*N1*( (a-X)*(b3->coordinateAux[2] - (1-position)*b1->coordinateAux[2] - position*b2->coordinateAux[2]) -c*(1-position)*(b2->coordinateAux[2] - b1->coordinateAux[2])) );
    
    //
    
    b4->forceAux[0] += -kDihed*DY*(b4->coordinateAux[0] - b3->coordinateAux[0]);
    
    b4->forceAux[1] += -kDihed*DY*(b4->coordinateAux[1] - b3->coordinateAux[1]);
    
    b4->forceAux[2] += -kDihed*DY*(b4->coordinateAux[2] - b3->coordinateAux[2]);
    
    
}