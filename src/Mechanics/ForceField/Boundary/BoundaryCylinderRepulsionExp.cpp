
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "BoundaryCylinderRepulsionExp.h"
#include "BoundaryCylinderRepulsion.h"

#include "BoundaryElement.h"
#include "Bead.h"

double BoundaryCylinderRepulsionExp::energy(double *coord, double *f, int *beadSet,
                                            double *krep, double *slen, int *nneighbors) {
    
    int nb, nc;
    double *coord1, R, r, U_i;
    double U = 0;
    auto Cumnc=0;
    auto beList = BoundaryElement::getBoundaryElements();
    nb = beList.size();
    
    for (int ib = 0; ib < nb; ib++) {
        
        auto be = beList[ib];
        nc = nneighbors[ib];
        
        for(int ic = 0; ic < nc; ic++) {
            
            coord1 = &coord[3 * beadSet[Cumnc + ic]];
            r = be->distance(coord1);
            
            R = -r / slen[Cumnc + ic];
            U_i = krep[Cumnc + ic] * exp(R);
            
//            std::cout<<r<<" "<<U_i<<endl;
            if(fabs(U_i) == numeric_limits<double>::infinity()
               || U_i != U_i || U_i < -1.0) {
                
                //set culprit and return
                BoundaryInteractions::_boundaryElementCulprit = be;
                ///TODO
                //BoundaryInteractions::_otherCulprit;
                
                return -1;
            }
            U += U_i;
        }
        Cumnc+=nc;
    }
    return U;
}

double BoundaryCylinderRepulsionExp::energy(double *coord, double *f, int *beadSet,
                                            double *krep, double *slen, int *nneighbors, double d) {
    
    int nb, nc;
    double *coord1, *force1, R, r, U_i;
    double U = 0;
    long Cumnc=0;
    auto beList = BoundaryElement::getBoundaryElements();
    nb = beList.size();
    
    for (int ib = 0; ib < nb; ib++) {
        
        auto be = beList[ib];
        nc = nneighbors[ib];
        
        for(int ic = 0; ic < nc; ic++) {
            
            coord1 = &coord[3 * beadSet[Cumnc + ic]];
            force1 = &f[3 * beadSet[Cumnc + ic]];
            
            r = be->stretchedDistance(coord1, force1, d);
            
            R = -r / slen[Cumnc + ic];
            
//            std::cout<<r<<" "<<krep[Cumnc+ic]<<endl;
            U_i = krep[Cumnc + ic] * exp(R);
            
            
            if(fabs(U_i) == numeric_limits<double>::infinity()
               || U_i != U_i || U_i < -1.0) {
                
                //set culprit and return
                BoundaryInteractions::_boundaryElementCulprit = be;
                ///TODO
                //BoundaryInteractions::_otherCulprit;
                
                return -1;
            }
            U += U_i;
        }
        Cumnc+=nc;
    }
    return U;
}



void BoundaryCylinderRepulsionExp::forces(double *coord, double *f, int *beadSet,
                                          double *krep, double *slen, int *nneighbors) {
    int nb, nc;
    double *coord1, *force1, R, r, f0;
    
    auto beList = BoundaryElement::getBoundaryElements();
    nb = beList.size();
    auto Cumnc=0;
    
    for (int ib = 0; ib < nb; ib++) {
        
        auto be = beList[ib];
        nc = nneighbors[ib];
//        std::cout<<"new "<<nb<<" "<<nc<<endl;
        for(int ic = 0; ic < nc; ic++) {
            
            coord1 = &coord[3 * beadSet[ Cumnc + ic]];
            force1 = &f[3 * beadSet[ Cumnc + ic]];
//            std::cout<<beadSet[0]<<" "<<beadSet[Cumnc + ic]<<endl;
            r = be->distance(coord1);
            auto norm = be->normal(coord1);
            
            R = -r / slen[Cumnc + ic];
            f0 = krep[Cumnc + ic] * exp(R);
            
//            std::cout<<"new "<<beadSet[ Cumnc + ic]<<" "<<force1[0]<<" "<<force1[1]<<" "<<force1[2]<<" "<<slen[Cumnc+ic]<<" "
//            <<krep[Cumnc+ic]<<" "<<f0<<endl;
            
            force1[0] += f0 *norm[0];
            force1[1] += f0 *norm[1];
            force1[2] += f0 *norm[2];
//            std::cout<<"BOUNDARY REPULSION EXP "<<force1[0]<<" "<<force1[1]<<" "<<force1[2]<<endl;
        }
        Cumnc+=nc;
    }
}

double BoundaryCylinderRepulsionExp::loadForces(double r, double kRep, double screenLength) {
    
    double R = -r/screenLength;
    return kRep * exp(R)/screenLength;

}
