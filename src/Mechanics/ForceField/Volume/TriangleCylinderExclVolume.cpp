
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

#include "TriangleCylinderExclVolume.h"

#include "MathFunctions.h"
using namespace mathfunc;

#include "TriangleCylinderBeadExclVolRepulsion.h"

#include "Triangle.h"
#include "Vertex.h"
#include "Cylinder.h"
#include "Bead.h"

template <class TriangleCylinderExclVolumeInteractionType>
double TriangleCylinderExclVolume<TriangleCylinderExclVolumeInteractionType>::computeEnergy(double d) {
    
    double U = 0;
    double U_i;

    for(auto t: Triangle::getTriangles()) {
        
        for(auto &c : _neighborList->getNeighbors(t)) {

            double kExVol = t->getMTriangle()->getExVolConst();

            // Use only 1st bead unless the cylinder is at plus end
            size_t numBeads = (c->isPlusEnd()? 2: 1);

            for(size_t idx = 0; idx < numBeads; ++idx) {
                Bead* b = (idx? c->getSecondBead(): c->getFirstBead());
                
                if (d == 0.0)
                    U_i = _FFType.energy(t, b, kExVol);
                else
                    U_i = _FFType.energy(t, b, kExVol, d);
                
                if(fabs(U_i) == numeric_limits<double>::infinity()
                || U_i != U_i || U_i < -1.0) {
                    
                    //set culprits and exit
                    _triangleCulprit = t;
                    _cylinderCulprit = c;
                    
                    return -1;
                }
                else
                    U += U_i;
            }
        }
    }
    
    return U;
}

template <class TriangleCylinderExclVolumeInteractionType>
void TriangleCylinderExclVolume<TriangleCylinderExclVolumeInteractionType>::computeForces() {

    for(auto t: Triangle::getTriangles()) {
        
        for(auto &c: _neighborList->getNeighbors(t)) {

            double kExVol = t->getMTriangle()->getExVolConst();

            // Use only 1st bead unless the cylinder is at plus end
            size_t numBeads = (c->isPlusEnd()? 2: 1);

            for(size_t idx = 0; idx < numBeads; ++idx) {
                Bead* b = (idx? c->getSecondBead(): c->getFirstBead());
            
                _FFType.forces(t, b, kExVol);
            }
        }
    }
}


template <class TriangleCylinderExclVolumeInteractionType>
void TriangleCylinderExclVolume<TriangleCylinderExclVolumeInteractionType>::computeForcesAux() {

    for(auto t: Triangle::getTriangles()) {
        
        for(auto &c: _neighborList->getNeighbors(t)) {

            double kExVol = t->getMTriangle()->getExVolConst();

            // Use only 1st bead unless the cylinder is at plus end
            size_t numBeads = (c->isPlusEnd()? 2: 1);

            for(size_t idx = 0; idx < numBeads; ++idx) {
                Bead* b = (idx? c->getSecondBead(): c->getFirstBead());
            
                _FFType.forcesAux(t, b, kExVol);
            }
        }
    }
}

template <class TriangleCylinderExclVolumeInteractionType>
void TriangleCylinderExclVolume<TriangleCylinderExclVolumeInteractionType>::computeLoadForces() {

    for(auto t: Triangle::getTriangles()) {
        
        for(auto &c: _neighborList->getNeighbors(t)) {

            double kExVol = t->getMTriangle()->getExVolConst();

            // potential acts on second cylinder bead if it is a plus  end
            // potential acts on first  cylinder bead if it is a minus end
            Bead* bd;
            Bead* bo;
            if(c->isPlusEnd()) {
                
                bd = c->getSecondBead();
                bo = c->getFirstBead();
                
                ///this normal is in the direction of polymerization
                auto normal = vector2Array<double, 3>(twoPointDirection(bo->coordinate, bd->coordinate));
                
                //array of coordinate values to update
                auto monSize = SysParams::Geometry().monomerSize[bd->getType()];
                auto cylSize = SysParams::Geometry().cylinderNumMon[bd->getType()];
                
                bd->lfip = 0;
                for (int i = 0; i < cylSize; i++) {
                    
                    array<double, 3> newCoord {{
                        bd->coordinate[0] + i * normal[0] * monSize,
                        bd->coordinate[1] + i * normal[1] * monSize,
                        bd->coordinate[2] + i * normal[2] * monSize
                    }};
                    
                    array<double, 3> loadForce = _FFType.loadForces(t, newCoord, kExVol);
                    bd->loadForcesP[bd->lfip++] += -dotProduct(normal, loadForce);
                }
                //reset lfi
                bd->lfip = 0;
            }
            
            if(c->isMinusEnd()) {
                
                bd = c->getFirstBead();
                bo = c->getSecondBead();
                
                ///this normal is in the direction of polymerization
                auto normal = vector2Array<double, 3>(twoPointDirection(bo->coordinate, bd->coordinate));
                
                //array of coordinate values to update
                auto monSize = SysParams::Geometry().monomerSize[bd->getType()];
                auto cylSize = SysParams::Geometry().cylinderNumMon[bd->getType()];
                
                
                bd->lfim = 0;
                for (int i = 0; i < cylSize; i++) {
                    
                    array<double, 3> newCoord {{
                        bd->coordinate[0] + i * normal[0] * monSize,
                        bd->coordinate[1] + i * normal[1] * monSize,
                        bd->coordinate[2] + i * normal[2] * monSize
                    }};
                    
                    array<double, 3> loadForce = _FFType.loadForces(t, newCoord, kExVol);
                    bd->loadForcesM[bd->lfim++] += -dotProduct(normal, loadForce);
                }
                //reset lfi
                bd->lfim = 0;
            }
            
        }
    }

}

///Template specializations
template double TriangleCylinderExclVolume<TriangleCylinderBeadExclVolRepulsion>::computeEnergy(double d);
template void TriangleCylinderExclVolume<TriangleCylinderBeadExclVolRepulsion>::computeForces();
template void TriangleCylinderExclVolume<TriangleCylinderBeadExclVolRepulsion>::computeForcesAux();
template void TriangleCylinderExclVolume<TriangleCylinderBeadExclVolRepulsion>::computeLoadForces();

