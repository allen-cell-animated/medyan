
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "Compartment.h"

#include "util/math/CuboidSlicing.hpp"
#include "GController.h"
#include "MathFunctions.h"
using namespace mathfunc;
#include "Visitor.h"
#include "Parser.h"
//REMOVE LATER
#include "ChemNRMImpl.h"

#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Triangle.h"
#include "Filament.h"
#include "Cylinder.h"
#include "GController.h"
#include <stdint.h>

Compartment& Compartment::operator=(const Compartment &other) {

    _species.clear();
    _internal_reactions.clear();
    _diffusion_reactions.clear();
    other.cloneSpecies(this);
    other.cloneReactions(this);
    _diffusion_rates = other._diffusion_rates;

    return *this;
}

bool Compartment::apply_impl(SpeciesVisitor &v) {
    for(auto &s : _species.species()) {
        v.visit(s.get());
    }
    return true;
}

bool Compartment::apply_impl(ReactionVisitor &v) {
    for(auto &r : _internal_reactions.reactions()) {
        v.visit(r.get());
    }
    return true;
}

void Compartment::computeNonSlicedVolumeArea() {
    auto sizex = SysParams::Geometry().compartmentSizeX;
    auto sizey = SysParams::Geometry().compartmentSizeY;
    auto sizez = SysParams::Geometry().compartmentSizeZ;

    _partialArea = {{sizey * sizez, sizey * sizez, sizex * sizez, sizex * sizez, sizex * sizey, sizex * sizey}};

    _volumeFrac = 1.0;
}

void Compartment::computeSlicedVolumeArea(SliceMethod sliceMethod) {

    switch(sliceMethod) {
    case SliceMethod::Membrane:
        {
            // The calculation requires the
            //  - The position calculation of triangles
            //  - The area calculation of triangles
            //  - The unit normal vector of triangles
            const size_t numTriangle = _triangles.size();
            if(numTriangle) {
                double sumArea = 0.0;
                Vec3 sumNormal {};
                Vec3 sumPos {};
                for(Triangle* t: _triangles) {
                    const auto& mesh = t->getParent()->getMesh();
                    const size_t ti = t->getTopoIndex();
                    const auto area = mesh.getTriangleAttribute(ti).gTriangle.area;
                    const auto& unitNormal = mesh.getTriangleAttribute(ti).gTriangle.unitNormal;
                    sumNormal += unitNormal * area;
                    sumPos += t->coordinate * area;
                    sumArea += area;
                }
                normalize(sumNormal);
                sumPos /= sumArea;

                PlaneCuboidSlicingResult res = PlaneCuboidSlicer() (
                    sumPos, sumNormal,
                    {
                        _coords[0] - SysParams::Geometry().compartmentSizeX * 0.5,
                        _coords[1] - SysParams::Geometry().compartmentSizeY * 0.5,
                        _coords[2] - SysParams::Geometry().compartmentSizeZ * 0.5
                    },
                    {{
                        SysParams::Geometry().compartmentSizeX,
                        SysParams::Geometry().compartmentSizeY,
                        SysParams::Geometry().compartmentSizeZ
                    }}
                );

                _volumeFrac = res.volumeIn / GController::getCompartmentVolume();
                _partialArea = res.areaIn;
            }
        }
        break;

    case SliceMethod::CylinderBoundary:
        {

            // The calculation requires the
            //  - The position calculation of triangles
            //  - The area calculation of triangles
            //  - The unit normal vector of triangles
            // ASSUMPTIONS:
            //  - This compartment is a CUBE
            //get compartment sizes in X,Y and the radius of cylinder
            auto sizex = SysParams::Geometry().compartmentSizeX;
            auto sizey = SysParams::Geometry().compartmentSizeY;
            auto sizez = SysParams::Geometry().compartmentSizeZ;
            auto r = SysParams::Boundaries().diameter / 2; //radius

            //get geometry center of the compartment
            auto x = _coords[0];
            auto y = _coords[1];

            auto leftx = x - sizex / 2;
            auto rightx = x + sizex / 2;
            auto lowy = y - sizey / 2;
            auto upy = y + sizey / 2;

            float pleft, pright, plow, pup, lleft, lright, llow, lup, VolumeIn;
            vector<float> edge;
            // edge_index = intersection points at left = 1, at right = 2, at low = 4 and at up = 5 in 2D;
            vector<int> edge_index;

            //1. find intersection points at left or right edges
            //if at lower or upper phase
            if(y < r){
                pleft = r - sqrt(r * r - (leftx - r) * (leftx - r));
                pright = r - sqrt(r * r - (rightx - r) * (rightx - r));

                //if the intersection is not inside the compartment, use the full compartlent size
                if(pleft > upy || pleft < lowy) lleft = sizey;
                else{
                    lleft = upy - pleft;
                    edge.push_back(lleft);
                    edge_index.push_back(1);
                }

                if(pright > upy || pright < lowy) lright = sizey;
                else{
                    lright = upy - pright;
                    edge.push_back(lright);
                    edge_index.push_back(2);

                }
            }
            else if(y > r){
                pleft = r + sqrt(r * r - (leftx - r) * (leftx - r));
                pright = r + sqrt(r * r - (rightx - r) * (rightx - r));

                //if the intersection is not inside the compartment, use the full compartlent size
                if(pleft > upy || pleft < lowy) lleft = sizey;
                else{
                    lleft = pleft - lowy;
                    edge.push_back(lleft);
                    edge_index.push_back(1);
                }

                if(pright > upy || pright < lowy) lright = sizey;
                else{
                    lright = pright - lowy;
                    edge.push_back(lright);
                    edge_index.push_back(2);
                }
            }
            else {
                cout<<"Even number of compartments in X or Y direction is not yet supportted."<<endl;
            }

            //1. find intersection points at lower or upper edges
            //if at left or right phase
            if(x < r){
                plow = r - sqrt(r * r - (lowy - r) * (lowy - r));
                pup = r - sqrt(r * r - (upy - r) * (upy - r));

                //if the intersection is not inside the compartment, use the full compartlent size
                if(plow > rightx || plow < leftx) llow = sizex;
                else{
                    llow = rightx - plow;
                    edge.push_back(llow);
                    edge_index.push_back(4);
                }

                if(pup > rightx || pup < leftx) lup = sizex;
                else{
                    lup = rightx - pup;
                    edge.push_back(lup);
                    edge_index.push_back(5);
                }
            }
            else if (x > r){
                plow = r + sqrt(r * r - (lowy - r) * (lowy - r));
                pup = r + sqrt(r * r - (upy - r) * (upy - r));

                //if the intersection is not inside the compartment, use the full compartlent size
                if(plow > rightx || plow < leftx) llow = sizex;
                else{
                    llow = plow - leftx;
                    edge.push_back(llow);
                    edge_index.push_back(4);
                }

                if(pup > rightx || pup < leftx) lup = sizex;
                else{
                    lup = pup - leftx;
                    edge.push_back(lup);
                    edge_index.push_back(5);
                }
            }
            else{
                cout<<"Even number of compartments in X or Y direction is not yet supportted."<<endl;
            }

            _partialArea = {{lleft * sizez, lright * sizez, llow * sizez, lup * sizez, sizex * sizey, sizex * sizey}};

            if(!areEqual(sizex,sizey))
                cout << "Volume calculation requires X dimension and Y dimension to be the same." << endl;

            float totalVol = sizex * sizey * sizez;

            //there are either 2 intersection points or 0 intersection points
            if(edge.size() == 2 && edge_index.size() == 2){
                //case 1, trapezoid
                if(abs(edge_index[0] - edge_index[1]) == 1)
                    _volumeFrac = 0.5 * (edge[0] + edge[1]) * sizex * sizez / totalVol;
                else if(edge_index[0] - edge_index[1] == 0)
                    cout <<"Intersection points are at the same edge!" << endl;
                //case 2, trangle
                else{
                    if(x < r && y < r){
                        if(edge_index[0] == 2 || edge_index[1] == 2)
                            _volumeFrac = 0.5 * edge[0] * edge[1] * sizez / totalVol;
                        else
                            _volumeFrac = 1 - 0.5 * edge[0] * edge[1] * sizez / totalVol;
                    }
                    else if(x > r && y < r){
                        if(edge_index[0] == 1 || edge_index[1] == 1)
                            _volumeFrac = 0.5 * edge[0] * edge[1] * sizez / totalVol;
                        else
                            _volumeFrac = 1 - 0.5 * edge[0] * edge[1] * sizez / totalVol;
                    }
                    else if(x < r && y > r){
                        if(edge_index[0] == 2 || edge_index[1] == 2)
                            _volumeFrac = 0.5 * edge[0] * edge[1] * sizez /totalVol;
                        else
                            _volumeFrac = 1 - 0.5 * edge[0] * edge[1] * sizez / totalVol;
                    }
                    else if(x > r && y > r){
                        if(edge_index[0] == 1 || edge_index[1] == 1)
                            _volumeFrac = 0.5 * edge[0] * edge[1] * sizez / totalVol;
                        else
                            _volumeFrac = 1 - 0.5 * edge[0] * edge[1] * sizez / totalVol;
                    }
                }
            }
            //case 3, no intersections.
            else if(edge.size() == 0 && edge_index.size() == 0){
                _volumeFrac = sizex * sizey * sizez / totalVol;
            }
            //case 4, two intersections points are the two vertices
            else if(edge.size() == 4 && edge_index.size() == 4){
                _volumeFrac = 0.5;
            }
            //case 5, only one intersection point is a vertex
            else if(edge.size() == 3 && edge_index.size() == 3){
                float a1;
                for(int i=0; i < 3; i++){
                    if(!areEqual(edge[i], 0.0) && !areEqual(edge[i], sizex))
                        a1 = edge[i];
                }
                _volumeFrac = 0.5 * a1 * sizex * sizez / totalVol;
            }
            else{
                cout <<"There are "<< edge.size() <<" intersection points for this compartment:"<< endl;
                cout << "x = " << _coords[0] << ", y = " << _coords[1] << ", z = " << _coords[2] <<endl;
                cout << "Something goes wrong!" << endl;
            }
        }
        break;
    }

}

vector<ReactionBase*> Compartment::generateDiffusionReactions(Compartment* C, bool outwardOnly) {
    // The compartment C and "this" must be neighbors of each other, and
    // "this" must be an active compartment.

    vector<ReactionBase*> rxns;

    // cout << "This compartment: x = " << _coords[0] << ", y = " << _coords[1] << ", z = " << _coords[2] <<endl;

    for(auto &sp_this : _species.species()) {
        int molecule = sp_this->getMolecule();
        float diff_rate = _diffusion_rates[molecule];
        if(diff_rate<0)  continue;

        if(C->isActivated()) {
            // Scale the diffusion rate according to the contacting areas
            size_t idxFwd = _neighborIndex.at(C), idxBwd = C->_neighborIndex.at(this);
            double scaleFactor = 0.5 * (_partialArea[idxFwd] + C->_partialArea[idxBwd]) / GController::getCompartmentArea()[idxFwd / 2];
            //double scaleFactor = 1.0;
            // cout << "To neighbor: x = " << C->_coords[0] << ", y = " << C->_coords[1] << ", z = " << C->_coords[2] <<endl;
            // cout << "scaleFactor = " << scaleFactor << endl;

            float actualDiffRate = diff_rate * scaleFactor;
            float volumeFrac = getVolumeFrac();
            // cout << "VolumeFraction = " << volumeFrac << endl;

            Species *sp_neighbour = C->_species.findSpeciesByMolecule(molecule);

            ReactionBase *R = new DiffusionReaction({sp_this.get(),sp_neighbour}, actualDiffRate, false, volumeFrac);
            this->addDiffusionReaction(R);
            rxns.push_back(R);

            if(!outwardOnly) {
                // Generate inward diffusion reaction
                ReactionBase* R = new DiffusionReaction({sp_neighbour, sp_this.get()}, actualDiffRate, false, C->getVolumeFrac());
                C->addDiffusionReaction(R);
                rxns.push_back(R);
            }
        }
    }

    return vector<ReactionBase*>(rxns.begin(), rxns.end());
}

vector<ReactionBase*> Compartment::generateAllDiffusionReactions(bool outwardOnly) {
    
    vector<ReactionBase*> rxns;

    if(_activated) {
        for (auto &C: _neighbours) {
            auto newRxns = generateDiffusionReactions(C, outwardOnly);
            rxns.insert(rxns.begin(), newRxns.begin(), newRxns.end());
        }
    }
    return vector<ReactionBase*>(rxns.begin(), rxns.end());
}

void Compartment::removeDiffusionReactions(ChemSim* chem, Compartment* C)
{
    //look for neighbor's diffusion reactions
    vector<ReactionBase*> to_remove;

    for(auto &r : C->_diffusion_reactions.reactions()) {

        auto rs = r.get()->rspecies()[1];
        if(rs->getSpecies().getParent() == this) {

            r->passivateReaction();

            chem->removeReaction(r.get());

            to_remove.push_back(r.get());
        }

    }

    //remove them
    for(auto &r : to_remove)
        C->_diffusion_reactions.removeReaction(r);

}

void Compartment::removeAllDiffusionReactions(ChemSim* chem) {

    //remove all diffusion reactions that this has ownership of
    for(auto &r : _diffusion_reactions.reactions()) {
        r->passivateReaction();
        chem->removeReaction(r.get());
    }

    _diffusion_reactions.clear();

    //remove neighboring diffusing reactions with this compartment
    for (auto &C: _neighbours)
        removeDiffusionReactions(chem, C);
}


void Compartment::transferSpecies(int i) {
    //i axis
    //0 X
    //1 Y
    //2 Z
    //3 all directions NOT IMPLEMENTED
    //get active neighbors
    vector<Compartment*> activeNeighbors;

    for(auto &neighbor : _neighbours){
        auto ncoord=neighbor->coordinates();

        if(neighbor->isActivated()){
            if(i==3) {
                activeNeighbors.push_back(neighbor);
                //Not implemented.
            }
            else if(mathfunc::twoPointDistance(ncoord,_coords)==(abs(_coords[i]-ncoord[i])))
                activeNeighbors.push_back(neighbor);
        }}

    assert(activeNeighbors.size() != 0
           && "Cannot transfer species to another compartment... no neighbors are active");
    if(i<3 && activeNeighbors.size()>1){
        cout<<"Error transferring species along an axis. More than 1 neighbor. Exiting. "<< endl;
        exit(EXIT_FAILURE);
    }

    //go through species
    Species* sp_neighbor;
    vector<Species*> sp_neighbors;

    for(auto &sp : _species.species()) {

        int copyNumber = sp->getN();
        auto nit = activeNeighbors.begin();

        if(sp->getFullName().find("Bound") == string::npos){
            while(copyNumber > 0) {
                sp->down();

                //choose a random active neighbor
                auto neighbor = *nit;

                sp_neighbor = neighbor->findSpeciesByName(sp->getName());

                //add to list if not already
                auto spit = find(sp_neighbors.begin(),
                                 sp_neighbors.end(),
                                 sp_neighbor);

                if(spit == sp_neighbors.end())
                    sp_neighbors.push_back(sp_neighbor);

                //increase copy number

                sp_neighbor->up();

                //reset if we've looped through
                if(++nit == activeNeighbors.end())
                    nit = activeNeighbors.begin();
                copyNumber--;

            }
        }

        //activate all reactions changed
        for(auto spn : sp_neighbors)
            spn->updateReactantPropensities();
        for(auto &sp : _species.species())
            sp->updateReactantPropensities();
    }
}

void Compartment::shareSpecies(int i) {
    //i axis
    //0 X
    //1 Y
    //2 Z
    //3 all directions
    //get active neighbors
    vector<Compartment*> activeNeighbors;

    for(auto &neighbor : _neighbours){
        auto ncoord=neighbor->coordinates();
    if(neighbor->isActivated()){
        if(i==3)
            activeNeighbors.push_back(neighbor);
        else if(mathfunc::twoPointDistance(ncoord,_coords)==(abs(_coords[i]-ncoord[i])))
        activeNeighbors.push_back(neighbor);
    }}

    assert(activeNeighbors.size() != 0
           && "Cannot share species to another compartment... no neighbors are active");
    if(i<3 && activeNeighbors.size()>1){
        cout<<"Error sharing species along an axis. More than 1 neighbor. Exiting."<< endl;
        exit(EXIT_FAILURE);
    }
    //go through species
    Species* sp_neighbor;
    vector<Species*> sp_neighbors;

    for(auto &sp : _species.species()) {
        auto nit = activeNeighbors.begin();
        auto neighbor = *nit;
        sp_neighbor = neighbor->findSpeciesByName(sp->getName());
        int copyNumber = sp_neighbor->getN();
        int lowerlimit = (int) sp_neighbor->getN()/2;
        if(sp->getFullName().find("Bound") == string::npos){
            while(copyNumber > lowerlimit) {
                sp_neighbor->down();

                //add to list if not already
                auto spit = find(sp_neighbors.begin(),
                                 sp_neighbors.end(),
                                 sp_neighbor);

                if(spit == sp_neighbors.end())
                    sp_neighbors.push_back(sp_neighbor);

                //increase copy number
                sp->up();
                //reset if we've looped through
                if(++nit == activeNeighbors.end())
                    nit = activeNeighbors.begin();
                neighbor = *nit;
                sp_neighbor = neighbor->findSpeciesByName(sp->getName());
                copyNumber--;

            }
        }

        //activate all reactions changed
        for(auto spn : sp_neighbors)
            spn->updateReactantPropensities();
        for(auto &sp : _species.species())
            sp->updateReactantPropensities();

    }
}

void Compartment::activate(ChemSim* chem) {
    /**************************************************************************
    The diffusion-reactions with the already activated neighbors would be added
    for both directions.
    **************************************************************************/
    
    assert(!_activated && "Compartment is already activated.");

    //set marker
    _activated = true;
    //add all diffusion reactions
    auto rxns = generateAllDiffusionReactions(false);

    for(auto &r : rxns) {
        chem->addReaction(r);
        r->activateReaction(); // Conditionally activate the new diffusion reactions
    }

    shareSpecies(SysParams::Boundaries().transfershareaxis); // FIXME: disable with membrane

}

void Compartment::updateActivation(ChemSim* chem) {
    double volumeFrac = getVolumeFrac();

    if(_activated) {
        // Update the reaction rates for diffusions in both directions
        for(auto& c: _neighbours) if(c->isActivated()) {
            // For any activated neighbor

            for(auto &sp_this : _species.species()) {
                int molecule = sp_this->getMolecule();
                float baseDiffRate = _diffusion_rates[molecule];
                if(baseDiffRate<0)  continue;

                Species *sp_neighbor = c->_species.findSpeciesByMolecule(molecule);

                // Scale the diffusion rate according to the contacting areas
                size_t idxFwd = _neighborIndex.at(c), idxBwd = c->_neighborIndex.at(this);
                double scaleFactor =
                    0.5 * (_partialArea[idxFwd] + c->_partialArea[idxBwd]) /
                    GController::getCompartmentArea()[idxFwd / 2];
                double actualDiffRate = baseDiffRate * scaleFactor;

                // Update outward reaction rate
                for(auto& r: _diffusion_reactions.reactions())
                    if(sp_this.get() == &r->rspecies()[0]->getSpecies() && sp_neighbor == &r->rspecies()[1]->getSpecies()) {
                        r->setVolumeFrac(volumeFrac);
                        r->setRateScaled(actualDiffRate);
                    }
                // We also update inward reaction rate here to ensure that neighbors are always on the same page.
                // Update inward reaction rate
                for(auto& r: c->_diffusion_reactions.reactions())
                    if(sp_this.get() == &r->rspecies()[1]->getSpecies() && sp_neighbor == &r->rspecies()[0]->getSpecies()) {
                        r->setVolumeFrac(c->getVolumeFrac());
                        r->setRateScaled(actualDiffRate);
                    }
            }

        }
    } else {
        activate(chem);
    }

    // Update the internal reaction rates
    for(auto& r: _internal_reactions.reactions()) {
        r->setVolumeFrac(volumeFrac);
        r->setRateScaled(r->getBareRate());
    }
}

void Compartment::deactivate(ChemSim* chem, bool init) {
    
    //assert no cylinders in this compartment
    assert((_cylinders.size() == 0)
           && "Compartment cannot be deactivated when containing active cylinders.");

    assert(_activated && "Compartment is already deactivated.");

    //set marker
    _activated = false;
    
    if(!init) transferSpecies(SysParams::Mechanics().transfershareaxis);
    removeAllDiffusionReactions(chem);
}

bool operator==(const Compartment& a, const Compartment& b) {
    if(a.numberOfSpecies()!=b.numberOfSpecies() or
       a.numberOfInternalReactions()!=b.numberOfInternalReactions())
        return false;

    if(typeid(a)!=typeid(b))
        return false;

    bool spec_bool = false;
    auto sit_pair = mismatch(a._species.species().begin(),
                             a._species.species().end(),
                             b._species.species().begin(),
            [](const unique_ptr<Species> &A, const unique_ptr<Species> &B)
            {return (*A)==(*B); });
    if(sit_pair.first==a._species.species().end())
        spec_bool=true;


    bool reac_bool = false;
    auto rit_pair = mismatch(a._internal_reactions.reactions().begin(),
                             a._internal_reactions.reactions().end(),
                             b._internal_reactions.reactions().begin(),
            [](const unique_ptr<ReactionBase> &A, const unique_ptr<ReactionBase> &B)
            {return (*A)==(*B);});
    if(rit_pair.first==a._internal_reactions.reactions().end())
        reac_bool=true;

    return spec_bool && reac_bool;
}
