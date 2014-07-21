//
//  MFilamentBending.cpp
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#include "MFilament.h"
#include "MBead.h"
#include "MCylinder.h"

using namespace std;
using namespace mathfunc;

double Filament::EnergyHarmonicBending(Cylinder* pc1, Cylinder* pc2 ){
    
    
    Bead* pb1 = pc1->GetFirstBead();
    Bead* pb2 = pc2->GetFirstBead();
    Bead* pb3 = pc2->GetSecondBead();
    

    
    double L1 = sqrt(ScalarProduct(pb1->coordinate, pb2->coordinate, pb1->coordinate, pb2->coordinate));
    double L2 = sqrt(ScalarProduct(pb2->coordinate, pb3->coordinate, pb2->coordinate, pb3->coordinate));
    
    double L1L2 = L1*L2;
    
    double l1l2 = ScalarProduct(pb1->coordinate, pb2->coordinate, pb2->coordinate, pb3->coordinate);
    
    double u = pc2->GetBendingConst() * ( 1 - l1l2 / L1L2 );
    
    pb1 = NULL;
    pb2 = NULL;
    pb3 = NULL;
    
    return u;
}

double Filament::EnergyHarmonicBending(Cylinder* pc1, Cylinder* pc2, double d ){
    
    Bead* pb1 = pc1->GetFirstBead();
    Bead* pb2 = pc2->GetFirstBead();
    Bead* pb3 = pc2->GetSecondBead();
    
    double L1 = sqrt(ScalarProductStretched(pb1->coordinate, pb1->force, pb2->coordinate, pb2->force, pb1->coordinate, pb1->force, pb2->coordinate, pb2->force, d));
    double L2 = sqrt(ScalarProductStretched(pb1->coordinate, pb1->force, pb2->coordinate, pb2->force, pb1->coordinate, pb1->force, pb2->coordinate, pb2->force, d));
    
    double L1L2 = L1*L2;
    
    double l1l2 = ScalarProductStretched(pb1->coordinate, pb1->force, pb2->coordinate, pb2->force, pb1->coordinate, pb1->force, pb2->coordinate, pb2->force, d);
    
    double u = pc2->GetBendingConst() * ( 1 - l1l2 / L1L2 );
    
    pb1 = NULL;
    pb2 = NULL;
    pb3 = NULL;
    
    return u;
}

void Filament::ForceHarmonicBending(Cylinder* pc1, Cylinder* pc2 ){
    
    Bead* pb1 = pc1->GetFirstBead();
    Bead* pb2 = pc2->GetFirstBead();
    Bead* pb3 = pc2->GetSecondBead();
    
    double k = pc2->GetBendingConst();
    
    double L1 = sqrt(ScalarProduct(pb1->coordinate, pb2->coordinate, pb1->coordinate, pb2->coordinate));
    double L2 = sqrt(ScalarProduct(pb2->coordinate, pb3->coordinate, pb2->coordinate, pb3->coordinate));
    double l1l2 = ScalarProduct(pb1->coordinate, pb2->coordinate, pb2->coordinate, pb3->coordinate);
    
    //invL = 1/L;
    double invL1 = 1/L1;
    double invL2 = 1/L2;
    double A = invL1*invL2;
    double B = l1l2*invL1*A*A*L2;
    double C = l1l2*invL2*A*A*L1;
    
    //force on i-1, f = k*(-A*l2 + B*l1):
    pb1->force[0] += k * ( (-pb3->coordinate[0] + pb2->coordinate[0])*A + (pb2->coordinate[0] - pb1->coordinate[0])*B );
    pb1->force[1] += k * ( (-pb3->coordinate[1] + pb2->coordinate[1])*A + (pb2->coordinate[1] - pb1->coordinate[1])*B );
    pb1->force[2] += k * ( (-pb3->coordinate[2] + pb2->coordinate[2])*A + (pb2->coordinate[2] - pb1->coordinate[2])*B );
    
    //force on i, f = k*(A*(l1-l2) - B*l1 + C*l2):
    pb2->force[0] += k *( (pb3->coordinate[0] - 2*pb2->coordinate[0] + pb1->coordinate[0])*A - (pb2->coordinate[0] - pb1->coordinate[0])*B + (pb3->coordinate[0] - pb2->coordinate[0])*C );
    pb2->force[1] += k *( (pb3->coordinate[1] - 2*pb2->coordinate[1] + pb1->coordinate[1])*A - (pb2->coordinate[1] - pb1->coordinate[1])*B + (pb3->coordinate[1] - pb2->coordinate[1])*C );
    pb2->force[2] += k *( (pb3->coordinate[2] - 2*pb2->coordinate[2] + pb1->coordinate[2])*A - (pb2->coordinate[2] - pb1->coordinate[2])*B + (pb3->coordinate[2] - pb2->coordinate[2])*C );
    
    //force on i-1, f = k*(A*l - B*l2):
    pb3->force[0] += k *( (pb2->coordinate[0] - pb1->coordinate[0])*A - (pb3->coordinate[0] - pb2->coordinate[0])*C );
    pb3->force[1] += k *( (pb2->coordinate[1] - pb1->coordinate[1])*A - (pb3->coordinate[1] - pb2->coordinate[1])*C );
    pb3->force[2] += k *( (pb2->coordinate[2] - pb1->coordinate[2])*A - (pb3->coordinate[2] - pb2->coordinate[2])*C );
    
    pb1 = NULL;
    pb2 = NULL;
    pb3 = NULL;
    
}

void Filament::ForceHarmonicBendingAux(Cylinder* pc1, Cylinder* pc2 ){
    
    
    Bead* pb1 = pc1->GetFirstBead();
    Bead* pb2 = pc2->GetFirstBead();
    Bead* pb3 = pc2->GetSecondBead();
    
    double k = pc2->GetBendingConst();
    
    double L1 = sqrt(ScalarProduct(pb1->coordinate, pb2->coordinate, pb1->coordinate, pb2->coordinate));
    double L2 = sqrt(ScalarProduct(pb2->coordinate, pb3->coordinate, pb2->coordinate, pb3->coordinate));
    double l1l2 = ScalarProduct(pb1->coordinate, pb2->coordinate, pb2->coordinate, pb3->coordinate);
    
    //invL = 1/L;
    double invL1 = 1/L1;
    double invL2 = 1/L2;
    double A = invL1*invL2;
    double B = l1l2*invL1*A*A*L2;
    double C = l1l2*invL2*A*A*L1;
    
    //force on i-1, f = k*(-A*l2 + B*l1):
    pb1->force[0] += k * ( (-pb3->coordinate[0] + pb2->coordinate[0])*A + (pb2->coordinate[0] - pb1->coordinate[0])*B );
    pb1->force[1] += k * ( (-pb3->coordinate[1] + pb2->coordinate[1])*A + (pb2->coordinate[1] - pb1->coordinate[1])*B );
    pb1->force[2] += k * ( (-pb3->coordinate[2] + pb2->coordinate[2])*A + (pb2->coordinate[2] - pb1->coordinate[2])*B );
    
    //force on i, f = k*(A*(l1-l2) - B*l1 + C*l2):
    pb2->forceAux[0] += k *( (pb3->coordinate[0] - 2*pb2->coordinate[0] + pb1->coordinate[0])*A - (pb2->coordinate[0] - pb1->coordinate[0])*B + (pb3->coordinate[0] - pb2->coordinate[0])*C );
    pb2->forceAux[1] += k *( (pb3->coordinate[1] - 2*pb2->coordinate[1] + pb1->coordinate[1])*A - (pb2->coordinate[1] - pb1->coordinate[1])*B + (pb3->coordinate[1] - pb2->coordinate[1])*C );
    pb2->forceAux[2] += k *( (pb3->coordinate[2] - 2*pb2->coordinate[2] + pb1->coordinate[2])*A - (pb2->coordinate[2] - pb1->coordinate[2])*B + (pb3->coordinate[2] - pb2->coordinate[2])*C );
    
    //force on i-1, f = k*(A*l - B*l2):
    pb3->forceAux[0] += k *( (pb2->coordinate[0] - pb1->coordinate[0])*A - (pb3->coordinate[0] - pb2->coordinate[0])*C );
    pb3->forceAux[1] += k *( (pb2->coordinate[1] - pb1->coordinate[1])*A - (pb3->coordinate[1] - pb2->coordinate[1])*C );
    pb3->forceAux[2] += k *( (pb2->coordinate[2] - pb1->coordinate[2])*A - (pb3->coordinate[2] - pb2->coordinate[2])*C );
    
    pb1 = NULL;
    pb2 = NULL;
    pb3 = NULL;
    
}