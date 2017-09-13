
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

#include "CylinderExclVolume.h"

#include "CylinderExclVolRepulsion.h"

#include "Cylinder.h"
#include "Bead.h"

#include "MathFunctions.h"
#include "cross_check.h"

using namespace mathfunc;

template <class CVolumeInteractionType>
int CylinderExclVolume<CVolumeInteractionType>::numInteractions;

template <class CVolumeInteractionType>
void CylinderExclVolume<CVolumeInteractionType>::vectorize() {
    
    //count interactions
    int nint = 0;
    
    for(auto ci : Cylinder::getCylinders()) {
        if(!ci->isFullLength()) continue;
        //do not calculate exvol for a non full length cylinder
        for(auto &cn : _neighborList->getNeighbors(ci))
        {
            if(!cn->isFullLength()||
               cn->getBranchingCylinder() == ci) continue;
            nint++;
        }
    }
    
    numInteractions = nint;
//    std::cout<<"NINT1 "<<nint<<endl;
    beadSet = new int[n * nint];
    krep = new double[nint];
    
    
    int nc = Cylinder::getCylinders().size();
    int i = 0;
    int Cumnc=0;
    for (i = 0; i < nc; i++) {
        auto ci = Cylinder::getCylinders()[i];
        if(!ci->isFullLength()) continue;

        int nn = _neighborList->getNeighbors(ci).size();
//        std::cout<<"Cylinder "<<i<<" "<<nn<<endl;
        for (int ni = 0; ni < nn; ni++) {
            
            auto cin = _neighborList->getNeighbors(ci)[ni];
            if(!cin->isFullLength()||
               cin->getBranchingCylinder() == ci) continue;
            
//            std::cout<<"N "<<n * (Cumnc + ni)<<endl;
//            std::cout<<beadSet[n * (Cumnc + ni)]<<endl;
//            std::cout<<ci->getFirstBead()->_dbIndex<<endl;
            beadSet[n * (Cumnc)] = ci->getFirstBead()->_dbIndex;
            beadSet[n * (Cumnc) + 1] = ci->getSecondBead()->_dbIndex;
            beadSet[n * (Cumnc) + 2] = cin->getFirstBead()->_dbIndex;
            beadSet[n * (Cumnc) + 3] = cin->getSecondBead()->_dbIndex;
//            std::cout<<Cumnc<<" "<<ni<<" "<<n * (Cumnc + ni)<<" "<<beadSet[n * (Cumnc + ni)]<<" "<<beadSet[n * (Cumnc + ni)+1]<<" "<<beadSet[n * (Cumnc + ni)+2]<<" "<<beadSet[n * (Cumnc + ni)+3]<<endl;
//            std::cout<<(Cumnc + ni)<<" "<<(Cumnc + ni)+1<<" "<<(Cumnc + ni)+2<<" "<<(Cumnc + ni)+3<<endl;
            krep[Cumnc] = ci->getMCylinder()->getExVolConst();
            Cumnc++;
        }
    }
}


template <class CVolumeInteractionType>
void CylinderExclVolume<CVolumeInteractionType>::deallocate() {
    
    delete beadSet;
    delete krep;
}


template <class CVolumeInteractionType>
double CylinderExclVolume<CVolumeInteractionType>::computeEnergy(double *coord, double *f, double d) {
    
    
    double U_i = 0;
    
    if (d == 0.0)
        U_i = _FFType.energy(coord, f, beadSet, krep);
    else
        U_i = _FFType.energy(coord, f, beadSet, krep, d);
//    std::cout<<"======================"<<endl;
#ifdef CROSSCHECK
    double U2 = 0;
    double U_ii;
    double nint=0;
    for(auto ci : Cylinder::getCylinders()) {
        
        //do not calculate exvol for a non full length cylinder
        if(!ci->isFullLength()) continue;
        
        for(auto &cn : _neighborList->getNeighbors(ci)) {
            
            //do not calculate exvol for a branching cylinder
            if(!cn->isFullLength() ||
               cn->getBranchingCylinder() == ci) continue;
            nint++;
            Bead* b1 = ci->getFirstBead();
            Bead* b2 = ci->getSecondBead();
            Bead* b3 = cn->getFirstBead();
            Bead* b4 = cn->getSecondBead();
//            std::cout<<b1->_dbIndex<<" "<<b2->_dbIndex<<" "<<b3->_dbIndex<<" "<<b4->_dbIndex<<endl;
            double kRepuls = ci->getMCylinder()->getExVolConst();
            
            if (d == 0.0)
                U_ii = _FFType.energy(b1, b2, b3, b4, kRepuls);
            else
                U_ii = _FFType.energy(b1, b2, b3, b4, kRepuls, d);
            
//            std::cout<<U_ii<<endl;
            if(fabs(U_i) == numeric_limits<double>::infinity()
               || U_ii != U_ii || U_ii < -1.0) {
                
                U2=-1;
            }
            else
                U2 += U_ii;
        }
    }
//    std::cout<<"NINT2 "<<nint<<endl;
    if(abs(U_i-U2)<=U2/100000000000)
        std::cout<<"E V YES "<<endl;
    else
    {   std::cout<<U_i<<" "<<U2<<endl;
        exit(EXIT_FAILURE);
    }

#endif
    return U_i;
}

template <class CVolumeInteractionType>
void CylinderExclVolume<CVolumeInteractionType>::computeForces(double *coord, double *f) {

    _FFType.forces(coord, f, beadSet, krep);
//    std::cout<<"==================="<<endl;
#ifdef CROSSCHECK
    double nint=0;
    for(auto ci : Cylinder::getCylinders()) {
        
        //do not calculate exvol for a non full length cylinder
        if(!ci->isFullLength()) continue;
        
        for(auto &cn : _neighborList->getNeighbors(ci)) {
            
            //do not calculate exvol for a branching cylinder
            if(!cn->isFullLength() ||
               cn->getBranchingCylinder() == ci) continue;
            nint++;

            Bead* b1 = ci->getFirstBead();
            Bead* b2 = ci->getSecondBead();
            Bead* b3 = cn->getFirstBead();
            Bead* b4 = cn->getSecondBead();
            std::cout<<b1->_dbIndex<<" "<<b2->_dbIndex<<" "<<b3->_dbIndex<<" "<<b4->_dbIndex<<endl;
            double kRepuls = ci->getMCylinder()->getExVolConst();
            
            if(cross_checkclass::Aux)
                _FFType.forcesAux(b1, b2, b3, b4, kRepuls);
            else
                _FFType.forces(b1, b2, b3, b4, kRepuls);

        }
    }
//    std::cout<<"NINT2 "<<nint<<endl;
        if(cross_checkclass::Aux){
            auto state=cross_check::crosscheckAuxforces(f);
            std::cout<<"F S+B+L+M+ +V YES "<<state<<endl;}
        else{
            auto state=cross_check::crosscheckforces(f);
            std::cout<<"F S+B+L+M+ +V YES "<<state<<endl;}
#endif
}

///Template specializations
template double CylinderExclVolume<CylinderExclVolRepulsion>::computeEnergy(double *coord, double *f, double d);
template void CylinderExclVolume<CylinderExclVolRepulsion>::computeForces(double *coord, double *f);
template void CylinderExclVolume<CylinderExclVolRepulsion>::vectorize();
template void CylinderExclVolume<CylinderExclVolRepulsion>::deallocate();


