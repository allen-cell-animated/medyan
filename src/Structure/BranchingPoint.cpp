
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

#include "BranchingPoint.h"

#include "SubSystem.h"
#include "Bead.h"
#include "Cylinder.h"
#include "Filament.h"
#include "ChemRNode.h"
#include "CompartmentGrid.h"

#include "GController.h"
#include "SysParams.h"
#include "MathFunctions.h"
#include "Rand.h"

using namespace mathfunc;

void BranchingPoint::updateCoordinate() {
    
    coordinate = midPointCoordinate(_c1->getFirstBead()->coordinate,
                                    _c1->getSecondBead()->coordinate,
                                    _position);
}

BranchingPoint::BranchingPoint(Cylinder* c1, Cylinder* c2,
                               short branchType, floatingpoint position)

    : Trackable(true,true), _c1(c1), _c2(c2), _position(position),
      _branchType(branchType), _branchID(_branchingPoints.getID()), _birthTime(tau()) {
    
    //Find compartment
    updateCoordinate();
        
    try {_compartment = GController::getCompartment(coordinate);}
    catch (exception& e) {
        cout << e.what();
        
        printSelf();
        
        exit(EXIT_FAILURE);
    }
        
    int pos = int(position * SysParams::Geometry().cylinderNumMon[c1->getType()]);
#ifdef CHEMISTRY
    _cBranchingPoint = unique_ptr<CBranchingPoint>(
    new CBranchingPoint(branchType, _compartment, c1->getCCylinder(), c2->getCCylinder(), pos));
    _cBranchingPoint->setBranchingPoint(this);
#endif
    
#ifdef MECHANICS
    _mBranchingPoint = unique_ptr<MBranchingPoint>(new MBranchingPoint(branchType));
    _mBranchingPoint->setBranchingPoint(this);
#endif
        
    //set the branching cylinder
    _c1->setBranchingCylinder(_c2);
}

BranchingPoint::~BranchingPoint() noexcept {
    
#ifdef MECHANICS
    //offset the branching cylinder's bead by a little for safety
    auto msize = SysParams::Geometry().monomerSize[_c1->getType()];
    
    vector<floatingpoint> offsetCoord =
    {(Rand::randInteger(0,1) ? -1 : +1) * Rand::randfloatingpoint(msize, 2 * msize),
     (Rand::randInteger(0,1) ? -1 : +1) * Rand::randfloatingpoint(msize, 2 * msize),
     (Rand::randInteger(0,1) ? -1 : +1) * Rand::randfloatingpoint(msize, 2 * msize)};
    
    auto b = _c2->getFirstBead();
    
    b->coordinate[0] += offsetCoord[0];
    b->coordinate[1] += offsetCoord[1];
    b->coordinate[2] += offsetCoord[2];
#endif
    
    
#ifdef CHEMISTRY
    //mark the correct species on the minus end of the branched
    //filament. If this is a filament species, change it to its
    //corresponding minus end. If a plus end, release a diffusing
    //or bulk species, depending on the initial reaction.
    CMonomer* m = _c2->getCCylinder()->getCMonomer(0);
    short speciesFilament = m->activeSpeciesFilament();
    
    //there is a filament species, mark its corresponding minus end
    if(speciesFilament != -1) {
        m->speciesMinusEnd(speciesFilament)->up();
        
        //unmark the filament and bound species
        m->speciesFilament(speciesFilament)->down();
        m->speciesBound(SysParams::Chemistry().brancherBoundIndex[_c1->getType()])->down();
    }
    //mark the free species instead
    else {
        //find the free species
        Species* speciesFilament = m->speciesFilament(m->activeSpeciesPlusEnd());
        
        string speciesName = SpeciesNamesDB::removeUniqueFilName(speciesFilament->getName());
        string speciesFirstChar = speciesName.substr(0,1);
        
        //find the free monomer, either bulk or diffusing
        Species* freeMonomer = nullptr;
        auto grid = _subSystem->getCompartmentGrid();
        
        Species* dMonomer  = _compartment->findSpeciesByName(speciesName);
        Species* dfMonomer = _compartment->findSpeciesByName(speciesFirstChar);
        
        Species* bMonomer  = grid->findSpeciesBulkByName(speciesName);
        Species* bfMonomer = grid->findSpeciesBulkByName(speciesFirstChar);
        
        //try diffusing
        if(dMonomer != nullptr) freeMonomer = dMonomer;
        // try bulk
        else if(bMonomer  != nullptr) freeMonomer = bMonomer;
        //diffusing, remove all but first char
        else if(dfMonomer != nullptr) freeMonomer = dfMonomer;
        //bulk, remove all but first char
        else if(bfMonomer != nullptr) freeMonomer = bfMonomer;
        //could not find. exit ungracefully
        else {
            cout << "In unbranching reaction, could not find corresponding " <<
                    "diffusing species of filament species " << speciesName <<
                    ". Exiting." << endl;
            exit(EXIT_FAILURE);
        }
            
        //remove the filament from the system
        Filament *bf = (Filament*)(_c2->getParent());
        _subSystem->removeTrackable<Filament>(bf);
        
        delete bf;
            
        //mark species, update reactions
        freeMonomer->up();
        freeMonomer->updateReactantPropensities();
    }
#endif
    //reset branching cylinder
    _c1->setBranchingCylinder(nullptr);
}

void BranchingPoint::updatePosition() {
    
#ifdef CHEMISTRY
    //update ccylinders
    _cBranchingPoint->setFirstCCylinder(_c1->getCCylinder());
    _cBranchingPoint->setSecondCCylinder(_c2->getCCylinder());
    
#endif
    //Find compartment
    updateCoordinate();
    
    Compartment* c;
    
    try {c = GController::getCompartment(coordinate);}
    catch (exception& e) {
        cout << e.what();
        
        printSelf();
        
        exit(EXIT_FAILURE);
    }
    
    if(c != _compartment) {
        _compartment = c;
#ifdef CHEMISTRY
        SpeciesBound* firstSpecies = _cBranchingPoint->getFirstSpecies();
        
        CBranchingPoint* clone = _cBranchingPoint->clone(c);
        setCBranchingPoint(clone);
        
        _cBranchingPoint->setFirstSpecies(firstSpecies);
#endif
    }
}
            
//Qin ----
void BranchingPoint::updateReactionRates() {
                
        //if no rate changers were defined, skip
        if(_unbindingChangers.empty()) return;
                
        //current force on branching point, use the total force
        floatingpoint fs = _mBranchingPoint->stretchForce;
        floatingpoint fb = _mBranchingPoint->bendingForce;
        floatingpoint fd = _mBranchingPoint->dihedralForce;
        floatingpoint ft = fs + fb + fd;
        floatingpoint force = max<floatingpoint>((floatingpoint)0.0, ft);
                
        //get the unbinding reaction
        ReactionBase* offRxn = _cBranchingPoint->getOffReaction();
                
        //change the rate
        float newRate = _unbindingChangers[_branchType]->changeRate(offRxn->getBareRate(), force);
        if(SysParams::RUNSTATE==false)
        {newRate=0.0;}
        offRxn->setRate(newRate);
        offRxn->updatePropensity();    
}
            
void BranchingPoint::printSelf() {
    
    cout << endl;
    
    cout << "BranchingPoint: ptr = " << this << endl;
    cout << "Branching type = " << _branchType << ", Branch ID = " << _branchID << endl;
    cout << "Coordinates = " << coordinate[0] << ", " << coordinate[1] << ", " << coordinate[2] << endl;
    
    cout << "Position on mother cylinder (floatingpoint) = " << _position << endl;
    cout << "Birth time = " << _birthTime << endl;
    
    cout << endl;
    
#ifdef CHEMISTRY
    cout << "Associated species = " << _cBranchingPoint->getFirstSpecies()->getName()
         << " , copy number = " << _cBranchingPoint->getFirstSpecies()->getN()
         << " , position on mother cylinder (int) = " << _cBranchingPoint->getFirstPosition() << endl;
#endif
    
    cout << endl;
    
    cout << "Associated cylinders (mother and branching): " << endl;
    _c1->printSelf();
    _c2->printSelf();
    
    cout << endl;
}
            
            

            
species_copy_t BranchingPoint::countSpecies(const string& name) {
    
    species_copy_t copyNum = 0;
    
    for(auto b : _branchingPoints.getElements()) {
        
        auto s = b->getCBranchingPoint()->getFirstSpecies();
        string sname = SpeciesNamesDB::removeUniqueFilName(s->getName());
        
        if(sname == name)
            copyNum += s->getN();
    }
    return copyNum;
}
            
vector<BranchRateChanger*> BranchingPoint::_unbindingChangers;

Database<BranchingPoint*> BranchingPoint::_branchingPoints;
            
