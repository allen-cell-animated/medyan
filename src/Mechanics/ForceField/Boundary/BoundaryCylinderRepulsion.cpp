
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

#include "BoundaryCylinderRepulsion.h"

#include "BoundaryCylinderRepulsionExp.h"
#include "BoundaryElement.h"

#include "Bead.h"
#include "Cylinder.h"

#include "MathFunctions.h"
#include "cross_check.h"

using namespace mathfunc;

template <class BRepulsionInteractionType>
void BoundaryCylinderRepulsion<BRepulsionInteractionType>::vectorize() {

    //count interactions
    int nint = 0;
    for (auto be: BoundaryElement::getBoundaryElements())
    {
        for(auto &c : _neighborList->getNeighbors(be))
        {
            if(c->isMinusEnd()) nint++;
            nint++;
        }
    }
//    std::cout<<"value of nint "<<nint<<endl;
    beadSet = new int[n * nint];
    krep = new double[nint];
    slen = new double[nint];
    
    
    int nbe = BoundaryElement::getBoundaryElements().size();
    int i = 0;
    int ni = 0;
    int bindex = 0;
    
    nneighbors = new int[nbe];
    auto cumnn=0;
    for (i = 0; i < nbe; i++) {

        auto be = BoundaryElement::getBoundaryElements()[i];
        auto nn = _neighborList->getNeighbors(be).size();
        
        nneighbors[i] = 0;
        auto idx=0;
        
        for (ni = 0; ni < nn; ni++) {
//            auto check=false;
            if(_neighborList->getNeighbors(be)[ni]->isMinusEnd())
            {
                bindex = _neighborList->getNeighbors(be)[ni]->getFirstBead()->_dbIndex;
                beadSet[cumnn+idx] = bindex;
                krep[cumnn+idx] = be->getRepulsionConst();
                slen[cumnn+idx] = be->getScreeningLength();
                idx++;
            }
                bindex = _neighborList->getNeighbors(be)[ni]->getSecondBead()->_dbIndex;
                beadSet[cumnn+idx] = bindex;
                krep[cumnn+idx] = be->getRepulsionConst();
                slen[cumnn+idx] = be->getScreeningLength();
                idx++;
            
            
//            if (_neighborList->getNeighbors(be)[ni]->isPlusEnd())
//            {bindex = _neighborList->getNeighbors(be)[ni]->getSecondBead()->_dbIndex;check=true;}
//            else if(_neighborList->getNeighbors(be)[ni]->isMinusEnd())
//            {bindex = _neighborList->getNeighbors(be)[ni]->getFirstBead()->_dbIndex;check=true;}
//                if(check){
//                    beadSet[cumnn+idx] = bindex;
//                    krep[cumnn+idx] = be->getRepulsionConst();
//                    slen[cumnn+idx] = be->getScreeningLength();
//                    idx++;
//                }
        }
        nneighbors[i]=idx;
        cumnn+=idx;
    }
}

template<class BRepulsionInteractionType>
void BoundaryCylinderRepulsion<BRepulsionInteractionType>::deallocate() {
    
    delete beadSet;
    delete krep;
    delete slen;
    delete nneighbors;
}

template <class BRepulsionInteractionType>
double BoundaryCylinderRepulsion<BRepulsionInteractionType>::computeEnergy(double *coord, double *f, double d) {

    
    double U_i;
    
    if (d == 0.0) {
        U_i = _FFType.energy(coord, f, beadSet, krep, slen, nneighbors);
    }
    else {
        U_i = _FFType.energy(coord, f, beadSet, krep, slen, nneighbors, d);
    }
//    std::cout<<"=================="<<endl;
#ifdef CROSSCHECK
    double U2 = 0;
    double U_ii;
    
    for (auto be: BoundaryElement::getBoundaryElements()) {
//        std::cout<<"neighbors "<<_neighborList->getNeighbors(be).size()<<endl;
        for(auto &c : _neighborList->getNeighbors(be)) {
            
            double kRep = be->getRepulsionConst();
            double screenLength = be->getScreeningLength();
            
            //potential acts on second bead unless this is a minus end
            Bead* bd;
            if(c->isMinusEnd()){
            bd = c->getFirstBead();
            
            if (d == 0.0)
            {U_ii =  _FFType.energy(
                                      bd, be->distance(bd->coordinate), kRep, screenLength);
//                std::cout<<be->distance(bd->coordinate)<<" "<<U_ii<<" ";
            }
            else
            {
                U_ii = _FFType.energy(
                                     bd, be->stretchedDistance(bd->coordinate, bd->force, d), kRep, screenLength);
//                std::cout<<be->stretchedDistance(bd->coordinate, bd->force, d)<<" "<<U_ii<<" ";
            }
            
            if(fabs(U_ii) == numeric_limits<double>::infinity()
               || U_ii != U_ii || U_ii < -1.0) {

                U2=-1;
                break;
            }
            else
                U2 += U_ii;
            }
            //---------------------------
            bd = c->getSecondBead();
            
            if (d == 0.0)
            {U_ii =  _FFType.energy(
                                    bd, be->distance(bd->coordinate), kRep, screenLength);
                //                std::cout<<be->distance(bd->coordinate)<<" "<<U_ii<<" ";
            }
            else
            {
                U_ii = _FFType.energy(
                                      bd, be->stretchedDistance(bd->coordinate, bd->force, d), kRep, screenLength);
                //                std::cout<<be->stretchedDistance(bd->coordinate, bd->force, d)<<" "<<U_ii<<" ";
            }
            
            
            if(fabs(U_ii) == numeric_limits<double>::infinity()
               || U_ii != U_ii || U_ii < -1.0) {
                
                U2=-1;
            }
            else
                U2 += U_ii;
            
        
        }
            }
    if(abs(U_i-U2)<=U2/100000000000)
        std::cout<<"E B YES "<<endl;
    else
    {   std::cout<<U_i<<" "<<U2<<endl;
        exit(EXIT_FAILURE);
    }

#endif
    return U_i;
}

template <class BRepulsionInteractionType>
void BoundaryCylinderRepulsion<BRepulsionInteractionType>::computeForces(double *coord, double *f) {
    for(auto b:Bead::getBeads())
        std::cout<<b->getID()<<" "<<b->_dbIndex<<endl;
    std::cout<<"====================="<<endl;
    for (auto be: BoundaryElement::getBoundaryElements()) {
        
        for(auto &c: _neighborList->getNeighbors(be)) {
            
            double kRep = be->getRepulsionConst();
            double screenLength = be->getScreeningLength();
            
            //potential acts on second cylinder bead unless this is a minus end
            if(c->isMinusEnd()) {
                std::cout<<"M "<<c->getFirstBead()->getID()<<" "<<c->getFirstBead()->_dbIndex<<" "<<
                c->getSecondBead()->_dbIndex<<" "<<c->getSecondBead()->getID()<<" "<<c->getMCylinder()->getLength()<<endl;
            }
            else if(c->isPlusEnd()) {
                std::cout<<"P "<<c->getFirstBead()->getID()<<" "<<c->getFirstBead()->_dbIndex<<" "<<
                c->getSecondBead()->_dbIndex<<" "<<c->getSecondBead()->getID()<<" "<<c->getMCylinder()->getLength()<<endl;
            }
            
        }
    }
    std::cout<<"++++++++++++++++++++++"<<endl;
    _FFType.forces(coord, f, beadSet, krep, slen, nneighbors);
    std::cout<<"====================="<<endl;
#ifdef CROSSCHECK
//    std::cout<<"old adjacency list"<<endl;
//    for (auto be: BoundaryElement::getBoundaryElements()) {
//        int count=0;
//        vector<int> b;
//        for(auto &c: _neighborList->getNeighbors(be)) {
//            if(c->isMinusEnd()) {count++; b.push_back(c->getFirstBead()->_dbIndex);}
//            else if(c->isPlusEnd()) {count++;b.push_back(c->getSecondBead()->_dbIndex);}
//            }
//        std::cout<<count<<endl;
//        for(auto i=0;i<count;i++)
//            std::cout<<b[i]<<" ";
//        std::cout<<endl;
//    }
//    for (auto be: BoundaryElement::getBoundaryElements()) {
//        for(auto &c: _neighborList->getNeighbors(be)) {
//            if(c->isMinusEnd()) {std::cout<<"old way "<<_neighborList->getNeighbors(be).size()<<" "<<c->getFirstBead()->_dbIndex<<" ";}
//            else if(c->isPlusEnd()) {std::cout<<"old way "<<_neighborList->getNeighbors(be).size()<<" "<<c->getFirstBead()->_dbIndex<<" ";}
//            std::cout<<endl;
//        }
//    }
    
    for (auto be: BoundaryElement::getBoundaryElements()) {
        
        for(auto &c: _neighborList->getNeighbors(be)) {

            double kRep = be->getRepulsionConst();
            double screenLength = be->getScreeningLength();
            
            //potential acts on second cylinder bead unless this is a minus end
            Bead* bd;
            if(c->isMinusEnd()) {
                bd = c->getFirstBead();
                auto normal = be->normal(bd->coordinate);
                if(cross_checkclass::Aux)
                    _FFType.forcesAux(bd, be->distance(bd->coordinate), normal, kRep, screenLength);
                else
                    _FFType.forces(bd, be->distance(bd->coordinate), normal, kRep, screenLength);
            }
            
                bd = c->getSecondBead();
                auto normal = be->normal(bd->coordinate);
                if(cross_checkclass::Aux)
                    _FFType.forcesAux(bd, be->distance(bd->coordinate), normal, kRep, screenLength);
                else
                    _FFType.forces(bd, be->distance(bd->coordinate), normal, kRep, screenLength);
        }
    }
    if(cross_checkclass::Aux)
        {    auto state=cross_check::crosscheckAuxforces(f);
            std::cout<<"F S+B+L+M+ +V+B YES "<<state<<endl;}
    else
    {    auto state=cross_check::crosscheckforces(f);
        std::cout<<"F S+B+L+M+ +V+B YES "<<state<<endl;}

    
#endif
}

template <class BRepulsionInteractionType>
void BoundaryCylinderRepulsion<BRepulsionInteractionType>::computeLoadForces() {
//    std::cout<<"BOUNDARY REPULSION DOES NOT USE VECTORIZED FORCES/COORDINATES"<<endl;
    for (auto be: BoundaryElement::getBoundaryElements()) {
        
        for(auto &c : _neighborList->getNeighbors(be)) {
            
            double kRep = be->getRepulsionConst();
            double screenLength = be->getScreeningLength();
            
            
            //potential acts on second cylinder bead unless this is a minus end
            Bead* bd;
            Bead* bo;
            if(c->isPlusEnd()) {
                
                bd = c->getSecondBead();
                bo = c->getFirstBead();
                
                ///this normal is in the direction of polymerization
                auto normal = normalizeVector(twoPointDirection(bo->coordinate, bd->coordinate));
                
                //array of coordinate values to update
                auto monSize = SysParams::Geometry().monomerSize[bd->getType()];
                auto cylSize = SysParams::Geometry().cylinderNumMon[bd->getType()];
                
                bd->lfip = 0;
                for (int i = 0; i < cylSize; i++) {
                    
                    auto newCoord = vector<double>{bd->coordinate[0] + i * normal[0] * monSize,
                        bd->coordinate[1] + i * normal[1] * monSize,
                        bd->coordinate[2] + i * normal[2] * monSize};
                    
                    double loadForce = _FFType.loadForces(be->distance(newCoord), kRep, screenLength);
                    bd->loadForcesP[bd->lfip++] += loadForce;
                }
                //reset lfi
                bd->lfip = 0;
            }
            
            if(c->isMinusEnd()) {
                
                bd = c->getFirstBead();
                bo = c->getSecondBead();
                
                ///this normal is in the direction of polymerization
                auto normal = normalizeVector(twoPointDirection(bo->coordinate, bd->coordinate));
                
                //array of coordinate values to update
                auto monSize = SysParams::Geometry().monomerSize[bd->getType()];
                auto cylSize = SysParams::Geometry().cylinderNumMon[bd->getType()];
                
                
                bd->lfim = 0;
                for (int i = 0; i < cylSize; i++) {
                    
                    auto newCoord = vector<double>{bd->coordinate[0] + i * normal[0] * monSize,
                        bd->coordinate[1] + i * normal[1] * monSize,
                        bd->coordinate[2] + i * normal[2] * monSize};
                    
                    double loadForce = _FFType.loadForces(be->distance(newCoord), kRep, screenLength);
                    bd->loadForcesM[bd->lfim++] += loadForce;
                }
                //reset lfi
                bd->lfim = 0;
            }
            
        }
        
    }
}

///Template specializations
template double BoundaryCylinderRepulsion<BoundaryCylinderRepulsionExp>::computeEnergy(double *coord, double *f, double d);
template void BoundaryCylinderRepulsion<BoundaryCylinderRepulsionExp>::computeForces(double *coord, double *f);
template void BoundaryCylinderRepulsion<BoundaryCylinderRepulsionExp>::computeLoadForces();
template void BoundaryCylinderRepulsion<BoundaryCylinderRepulsionExp>::vectorize();
template void BoundaryCylinderRepulsion<BoundaryCylinderRepulsionExp>::deallocate();

