
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
#include "MathFunctions.h"
#include "Visitor.h"
#include "Parser.h"
//REMOVE LATER
#include "ChemNRMImpl.h"
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

void Compartment::getNonSlicedVolumeArea() {
    auto sizex = SysParams::Geometry().compartmentSizeX;
    auto sizey = SysParams::Geometry().compartmentSizeY;
    auto sizez = SysParams::Geometry().compartmentSizeZ;

    _partialArea = {{sizey * sizez, sizey * sizez, sizex * sizez, sizex * sizez, sizex * sizey, sizex * sizey}};

    _partialVolume = 1.0;
}

//Calculates volume fraction
void Compartment::getSlicedVolumeArea() {
    // The calculation requires the
    //  - The position calculation of triangles
    //  - The area calculation of triangles
    //  - The unit normal vector of triangles
    // ASSUMPTIONS:
    //  - This compartment is a CUBE
//    size_t numTriangle = _triangles.size();
//    if(numTriangle) {
//        double sumArea = 0.0;
//        array<double, 3> sumNormal {};
//        array<double, 3> sumPos {};
//        for(Triangle* t: _triangles) {
//            double area = t->getGTriangle()->getArea();
//            vectorIncrease(sumNormal, vectorMultiply(t->getGTriangle()->getUnitNormal(), area));
//            vectorIncrease(sumPos, vectorMultiply(t->coordinate, area));
//            sumArea += area;
//        }
//        double oneOverSumArea = 1.0 / sumArea;
//        vectorExpand(sumNormal, oneOverSumArea);
//        vectorExpand(sumPos, oneOverSumArea);
//
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
            _partialVolume = 0.5 * (edge[0] + edge[1]) * sizex * sizez / totalVol;
        else if(edge_index[0] - edge_index[1] == 0)
            cout <<"Intersection points are at the same edge!" << endl;
        //case 2, trangle
        else{
            if(x < r && y < r){
                if(edge_index[0] == 2 || edge_index[1] == 2)
                    _partialVolume = 0.5 * edge[0] * edge[1] * sizez / totalVol;
                else
                    _partialVolume = 1 - 0.5 * edge[0] * edge[1] * sizez / totalVol;
            }
            else if(x > r && y < r){
                if(edge_index[0] == 1 || edge_index[1] == 1)
                    _partialVolume = 0.5 * edge[0] * edge[1] * sizez / totalVol;
                else
                    _partialVolume = 1 - 0.5 * edge[0] * edge[1] * sizez / totalVol;
            }
            else if(x < r && y > r){
                if(edge_index[0] == 2 || edge_index[1] == 2)
                    _partialVolume = 0.5 * edge[0] * edge[1] * sizez /totalVol;
                else
                    _partialVolume = 1 - 0.5 * edge[0] * edge[1] * sizez / totalVol;
            }
            else if(x > r && y > r){
                if(edge_index[0] == 1 || edge_index[1] == 1)
                    _partialVolume = 0.5 * edge[0] * edge[1] * sizez / totalVol;
                else
                    _partialVolume = 1 - 0.5 * edge[0] * edge[1] * sizez / totalVol;
            }
        }
    }
    //case 3, no intersections.
    else if(edge.size() == 0 && edge_index.size() == 0){
        _partialVolume = sizex * sizey * sizez / totalVol;
    }
    //case 4, two intersections points are the two vertices
    else if(edge.size() == 4 && edge_index.size() == 4){
        _partialVolume = 0.5;
    }
    //case 5, only one intersection point is a vertex
    else if(edge.size() == 3 && edge_index.size() == 3){
        float a1;
        for(int i=0; i < 3; i++){
            if(!areEqual(edge[i], 0.0) && !areEqual(edge[i], sizex))
                a1 = edge[i];
        }
        _partialVolume = 0.5 * a1 * sizex * sizez / totalVol;
    }
    else{
        cout <<"There are "<< edge.size() <<" intersection points for this compartment:"<< endl;
        cout << "x = " << _coords[0] << ", y = " << _coords[1] << ", z = " << _coords[2] <<endl;
        cout << "Something goes wrong!" << endl;
    }




//        _partialVolume = res.volumeIn;
//        _partialArea = res.areaIn;
//    }
}

//TODO
vector<ReactionBase*> Compartment::generateDiffusionReactions(Compartment* C) {
    // The compartment C and "this" must be neighbors of each other, and
    // "this" must be an active compartment.

    vector<ReactionBase*> rxns;

    cout << "This compartment: x = " << _coords[0] << ", y = " << _coords[1] << ", z = " << _coords[2] <<endl;

    for(auto &sp_this : _species.species()) {
        int molecule = sp_this->getMolecule();
        float diff_rate = _diffusion_rates[molecule];
        if(diff_rate<0)  continue;

        if(C->isActivated()) {
            // Scale the diffusion rate according to the contacting areas
            size_t idxFwd = _neighborIndex.at(C), idxBwd = C->_neighborIndex.at(this);
            double scaleFactor = 0.5 * (_partialArea[idxFwd] + C->_partialArea[idxBwd]) / GController::getCompartmentArea()[idxFwd / 2];
            //double scaleFactor = 1.0;
            cout << "To neighbor: x = " << C->_coords[0] << ", y = " << C->_coords[1] << ", z = " << C->_coords[2] <<endl;
            cout << "scaleFactor = " << scaleFactor << endl;

            float actualDiffRate = diff_rate * scaleFactor;
            float volumeFrac = getVolumeFrac();
            cout << "VolumeFraction = " << volumeFrac << endl;

            Species *sp_neighbour = C->_species.findSpeciesByMolecule(molecule);
            //Diffusion reaction from "this" compartment to C.
            ReactionBase *R = new DiffusionReaction({sp_this.get(),sp_neighbour}, actualDiffRate, false, volumeFrac);
            this->addDiffusionReaction(R);
            rxns.push_back(R);

        }
    }

    return vector<ReactionBase*>(rxns.begin(), rxns.end());
}

//Diffusion is now scaled directly in Compartment::generateDiffusionReactions(Compartment* C).
vector<ReactionBase*> Compartment::generateScaleDiffusionReactions(Compartment* C)
{
    vector<ReactionBase*> rxns;

    cout << "neighbor: x = " << C->_coords[0] << ", y = " << C->_coords[1] <<endl;
    auto factor = generateScaleFactor(C);
    cout << "factor = " << factor << endl;

    for(auto &sp_this : _species.species()) {
        int molecule = sp_this->getMolecule();
        float diff_rate = _diffusion_rates[molecule];
        if(diff_rate<0)  continue;

        if(C->isActivated()) {
            Species *sp_neighbour = C->_species.findSpeciesByMolecule(molecule);

            auto diff_rate_s = diff_rate * factor;

            ReactionBase *R = new DiffusionReaction({sp_this.get(),sp_neighbour},diff_rate_s);
            this->addDiffusionReaction(R);
            rxns.push_back(R);
        }


    }


    return vector<ReactionBase*>(rxns.begin(), rxns.end());
}

//Generate a scaling factor for diffusion constant. For cylinder with 1 compartment in Z direction only
double Compartment::generateScaleFactor(Compartment* C)
{
    vector<ReactionBase*> rxns;

    //get compartment sizes in X,Y and the radius of cylinder
    auto lx = SysParams::Geometry().compartmentSizeX;
    auto ly = SysParams::Geometry().compartmentSizeY;
    auto lz = SysParams::Geometry().compartmentSizeZ;
    auto r = SysParams::Boundaries().diameter / 2; //radius
    //float c1;

    if((_coords[0] - lx/2) < r && (_coords[0] + lx/2) > r) {
        cout << "Diffusion Scaling failed" << endl;
        return 1.0;
    }

    if((_coords[1] - ly/2) < r && (_coords[1] + ly/2) > r) {
        cout << "Diffusion Scaling failed" << endl;
        return 1.0;
    }

    //get geometry center of the compartment
    auto x = _coords[0];
    auto y = _coords[1];
    //get geometry center of the neighbor compartment
    auto nx = C->_coords[0];
    auto ny = C->_coords[1];
    //c1 is the intersection line between compartment
    double c1;
    //c2 is the intersection between boundary and c1
    double c2;

    //scale diffusion rate based on compartment area
    //1. find the location of the neighbor compartment
    //if transport along x axis
    if(ny == y) {

        //2. calculate the interection point
        //if at lower half of the system
        if(y < r) {
            //if transport to the left neighbor
            if(nx < x) c1 = x - lx/2;
            //if transport to the right neighbor
            else c1 = x + lx/2;
            c2 = r - sqrt(r * r - (c1 - r) * (c1 - r));

            //3. calculate scaling factor
            //check if intersection is within this compartment
            if(c2 < (y + ly/2) && c2 > (y - ly/2)) {
                double factor = (y + ly/2 - c2) / ly;
                return factor;
            }
            else return 1;

        }
        //if at upper part
        else {
            //at left
            if(nx < x) c1 = x - lx/2;

            else c1 = x + lx/2; //right

            c2 = r + sqrt(r * r - (c1 - r) * (c1 - r));


            //3. calculate scaling factor
            if(c2 < (y + ly/2) && c2 > (y - ly/2)) {
                float factor = (c2 - y + ly/2) / ly;
                return factor;
            }
            else return 1;

        }
    }

    else if(nx == x){
        //2. calculate the interection point
        //if at left part
        if(x < r) {
            //if at lower
            if(ny < y) c1 = y - ly/2;

            //if at upper
            else c1 = y + ly/2;

            c2 = r - sqrt(r * r - (c1 - r) * (c1 - r));

            //3. calculate scaling factor
            //check if interaction is within compartment
            if(c2 < (x + lx/2) && c2 > (x - lx/2)) {
                float factor = (_coords[0] + lx/2 - c2) / lx;
                return factor;
            }
            else return 1;

        }
        //if at right part
        else {
            //at lower
            if(ny < y) c1 = y - ly/2;

            else c1 = y + ly/2; //right

            c2 = r + sqrt(r * r - (c1 - r) * (c1 - r));

            //3. calculate scaling factor
            if(c2 < (x + lx/2) && c2 > (x - lx/2)) {
                float factor = (c2 - x + lx/2) / lx;
                return factor;
            }
            else return 1;

        }

    }

}


vector<ReactionBase*> Compartment::generateAllDiffusionReactions() {

    vector<ReactionBase*> rxns;

    if(_activated) {
        for (auto &C: _neighbours) {

            auto newRxns = generateDiffusionReactions(C);
            rxns.insert(rxns.begin(), newRxns.begin(), newRxns.end());
        }
    }
    return vector<ReactionBase*>(rxns.begin(), rxns.end());
}

vector<ReactionBase*> Compartment::generateAllpairsDiffusionReactions() {

    vector<ReactionBase*> rxns;

    if(_activated) {
        for (auto &C: _neighbours) {
            if(C->isActivated()){
                //generates diffusion reactions both from and to the chosen compartment
            auto newRxns = generateDiffusionReactions(C);
                 rxns.insert(rxns.begin(), newRxns.begin(), newRxns.end());
            newRxns = C->generateDiffusionReactions(this);
                 rxns.insert(rxns.begin(), newRxns.begin(), newRxns.end());

            }
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

    assert(!_activated && "Compartment is already activated.");

    //set marker
    _activated = true;
    //add all diffusion reactions
    auto rxns = generateAllpairsDiffusionReactions();
    for(auto &r : rxns) chem->addReaction(r);
    shareSpecies(SysParams::Boundaries().transfershareaxis);

    for (auto &C: _neighbours){
        if(C->isActivated()){
            for(auto &r : C->_diffusion_reactions.reactions()) {
                auto rs = r.get()->rspecies()[1];//product
                if(rs->getSpecies().getParent() == this) {
                    auto rs1 = r.get()->rspecies()[0];
                    if(rs1->getN()>0 && r->isPassivated()){
                        r->activateReaction();
                    }
        }
    }
        }
    }

}

void Compartment::deactivate(ChemSim* chem) {

    //assert no cylinders in this compartment
    assert((_cylinders.size() == 0)
           && "Compartment cannot be deactivated when containing active cylinders.");

    assert(_activated && "Compartment is already deactivated.");

    //set marker
    _activated = false;

    transferSpecies(SysParams::Boundaries().transfershareaxis);
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
