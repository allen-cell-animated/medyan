
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

#include "Structure/Cylinder.h"

#include "SubSystem.h"
#include "CController.h"
#include "ChemManager.h"
#include "ChemRNode.h"

#include "Filament.h"
#include "Bead.h"

#include "GController.h"
#include "MathFunctions.h"

using namespace mathfunc;

//static vars needed to vectorize on-the-fly
int Cylinder::maxcindex = 0;
int Cylinder::vectormaxsize = 0;
int Cylinder::Ncyl = 0;
//bool Cylinder::triggercylindervectorization = false;
vector<int> Cylinder::removedcindex;//vector of bead indices that were once alloted to other

// Renew structured data for all cylinders
// This function will not invalidate stable indices
void Cylinder::updateData() {
    for(auto c : getCylinders()) {
        const auto si = c->getStableIndex();
        auto& data = getDbData().value[si];

        data.filamentId = static_cast<Filament*>(c->getParent())->getId();
        data.positionOnFilament = c->getPosition();
        data.compartmentId = c->getCompartment()->getId();
        data.beadCoord[0] = c->getFirstBead()->coordinate();
        data.beadCoord[1] = c->getSecondBead()->coordinate();
        data.coord = vector2Vec<3, double>(c->coordinate);
        data.type = c->getType();
        data.id = c->getId();

#ifdef CHEMISTRY
        data.chemCylinder = c->getCCylinder();
#endif
    }
}

// beads but are free to be reallocated now.
void Cylinder::revectorize(cylinder* cylindervec){
    int i = 0;
    Ncyl = numCylinders() - 1; // New cylinder is stored BEFORE this function, so temporarily exclude last
    for(size_t idx = 0; idx < Ncyl; ++idx){
        auto cyl = getElements()[idx];
        //set _dcIndex
        cyl->_dcIndex = i; // Here _dcIndex == getIndex()
        //copy attributes to a structure array
        cylindervec[i].filamentID = dynamic_cast<Filament*>(cyl->getParent())->getId();
        cylindervec[i].filamentposition = cyl->getPosition();
        cylindervec[i].beads[0] = cyl->getFirstBead();
        cylindervec[i].beads[1] = cyl->getSecondBead();
        cylindervec[i].cmpID = cyl->getCompartment()->getId();
        cylindervec[i].cindex = i;
        auto coord = cyl->coordinate;
        cylindervec[i].coord[0] = coord[0];
        cylindervec[i].coord[1] = coord[1];
        cylindervec[i].coord[2] = coord[2];
        cylindervec[i].type = cyl->getType();
        cylindervec[i].ID = cyl->getId();
        i++;
    }
    removedcindex.clear();
    maxcindex = Ncyl;
}

void  Cylinder::copytoarrays() {
    long i =getStableIndex();
    cylinder* cylindervec = CUDAcommon::serlvars.cylindervec;
    //copy attributes to a structure array
    cylindervec[i].filamentID = dynamic_cast<Filament*>(this->getParent())->getId();
    cylindervec[i].filamentposition = _position;
    cylindervec[i].beads[0] = _b1;
    cylindervec[i].beads[1] = _b2;
    cylindervec[i].cmpID = _compartment->getId();
    cylindervec[i].cindex = i;
    cylindervec[i].type = _type;
    cylindervec[i].ID = getId();
    //update coordinate in updatecoordinate
/*    auto coord = coordinate;
    cylindervec[i].coord[0] = coord[0];
    cylindervec[i].coord[1] = coord[1];
    cylindervec[i].coord[2] = coord[2];*/

    //other arrays needed
/*    ccylindervec[i] = _cCylinder.get();
    cylinderpointervec[i] = this;*/
}

void Cylinder::resetarrays() {
    cylinder* cylindervec = CUDAcommon::serlvars.cylindervec;
    resetcylinderstruct(cylindervec, _dcIndex);
}

void Cylinder::updateCoordinate() {
    coordinate = midPointCoordinate(_b1->vcoordinate(), _b2->vcoordinate(), 0.5);
    //update the coordiante in cylinder structure.

    Cylinder::getDbData().value[getStableIndex()].coord = 0.5 * (_b1->coordinate() + _b2->coordinate());
}


Cylinder::Cylinder(Composite* parent, Bead* b1, Bead* b2, short type, int position,
                   bool extensionFront, bool extensionBack, bool initialization)

    : Trackable(true, true, true, false),
      _b1(b1), _b2(b2), _type(type), _position(position),
      db_type(CylinderInfoData::CylinderInfo {}) {
    
    parent->addChild(unique_ptr<Component>(this));
    //revectorize if needed
    revectorizeifneeded();

    Ncyl = getElements().size() - 1;
    //check if you need to revectorize.
    cylinder* cylindervec = CUDAcommon::serlvars.cylindervec;
    //copy attributes to a structure array
    cylindervec[_dcIndex].filamentID = dynamic_cast<Filament*>(this->getParent())->getId();
    cylindervec[_dcIndex].filamentposition = _position;
    cylindervec[_dcIndex].beads[0] = _b1;
    cylindervec[_dcIndex].beads[1] = _b2;

    //Set coordinate
    updateCoordinate();

    try {_compartment = GController::getCompartment(coordinate);}
    catch (exception& e) {
        cout << e.what() << endl;
        exit(EXIT_FAILURE);
    }

   //add to compartment
   _compartment->addCylinder(this);
          
#ifdef MECHANICS
          //set eqLength according to cylinder size
          
    double eqLength  = twoPointDistance(b1->vcoordinate(), b2->vcoordinate());
    if(!SysParams::RUNSTATE) //RESTARTPHASE
    {
        int nummonomers = (int) round(eqLength/ SysParams::Geometry().monomerSize[type]);
        double tpd = eqLength;
              
        if(nummonomers ==0){
            eqLength = SysParams::Geometry().monomerSize[type];
        }
        else{
            eqLength = (nummonomers) * SysParams::Geometry().monomerSize[type];
            double mindis = abs(tpd - eqLength);

            for(auto i=nummonomers-1;i<=min(nummonomers+1, SysParams::Geometry().cylinderNumMon[type]);i++){
                if(mindis > abs(tpd - i * SysParams::Geometry().monomerSize[type]))
                {
                    eqLength = i * SysParams::Geometry().monomerSize[type];
                    mindis = abs(tpd - eqLength);
                }
            }
        }
              
            
    }
    _mCylinder = unique_ptr<MCylinder>(new MCylinder(_type, eqLength));
    _mCylinder->setCylinder(this);
#endif
    
#ifdef CHEMISTRY
    _cCylinder = unique_ptr<CCylinder>(new CCylinder(_compartment, this));
    _cCylinder->setCylinder(this);

    //copy further components to the array
    cylindervec[_dcIndex].cmpID = _compartment->getId();
    cylindervec[_dcIndex].cindex = _dcIndex;
    cylindervec[_dcIndex].type = _type;
    cylindervec[_dcIndex].ID = getId();

    //init using chem manager
    _chemManager->initializeCCylinder(_cCylinder.get(), extensionFront,
                                      extensionBack, initialization);
#endif

    // Update the stored data
    auto& data = getDbData().value[getStableIndex()];

    data.filamentId = static_cast<Filament*>(getParent())->getId();
    data.positionOnFilament = _position;
    data.compartmentId = _compartment->getId();
    data.beadCoord[0] = _b1->coordinate();
    data.beadCoord[1] = _b2->coordinate();
    data.coord = vector2Vec<3, double>(coordinate);
    data.type = getType();
    data.id = getId();

#ifdef CHEMISTRY
    data.chemCylinder = getCCylinder();
#endif

}

Cylinder::~Cylinder() noexcept {
    //remove from compartment
    _compartment->removeCylinder(this);
    
}

/// Get filament type
int Cylinder::getType() {return _type;}

void Cylinder::updatePosition() {

    //check if were still in same compartment, set new position
    updateCoordinate();

    Compartment* c;
    try {c = GController::getCompartment(coordinate);}
    catch (exception& e) {
        cout << e.what();
        
        printSelf();
        
        exit(EXIT_FAILURE);
    }
    
    if(c != _compartment) { 

#ifdef CHEMISTRY
        auto oldCompartment = _compartment;
        auto newCompartment = c;
#endif
        
        //remove from old compartment, add to new
        _compartment->removeCylinder(this);
        _compartment = c;
        _compartment->addCylinder(this);

#ifdef CHEMISTRY
        auto oldCCylinder = _cCylinder.get();
        
        //Remove old ccylinder from binding managers
        for(auto &manager : oldCompartment->getFilamentBindingManagers()) {
#ifdef NLORIGINAL
            manager->removePossibleBindings(oldCCylinder);
#endif
#ifdef NLSTENCILLIST
            manager->removePossibleBindingsstencil(oldCCylinder);
#endif
        }

        //clone and set new ccylinder
        CCylinder* clone = _cCylinder->clone(c);
        setCCylinder(clone);
        
        auto newCCylinder = _cCylinder.get();

        //change both CCylinder and Compartment ID in the vector
        auto& data = getDbData().value[getStableIndex()];
        data.compartmentId = _compartment->getId();
        data.chemCylinder = newCCylinder;
        
        //Add new ccylinder to binding managers
        for(auto &manager : newCompartment->getFilamentBindingManagers()){
#ifdef NLORIGINAL
            manager->addPossibleBindings(newCCylinder);
#endif
#ifdef NLSTENCILLIST
            //This directs call to Hybrid Binding Manager.
            manager->addPossibleBindingsstencil(newCCylinder);
#endif
        }
    }
#endif
    
#ifdef MECHANICS
    //update length
    _mCylinder->setLength(twoPointDistance(_b1->vcoordinate(),
                                           _b2->vcoordinate()));
#endif

}

/// @note -  The function uses the bead load force to calculate this changed rate.
/// If there is no force on the beads the reaction rates are set to the bare.

void Cylinder::updateReactionRates() {
    
    double force;
    
    //if no rate changer was defined, skip
    if(_polyChanger.empty()) return;
    
    //load force from front (affects plus end polymerization)
    if(_plusEnd) {
        
        //get force of front bead
        force = _b2->getLoadForcesP();
        
        //change all plus end polymerization rates
        for(auto &r : _cCylinder->getInternalReactions()) {
            float newRate;
            if(r->getReactionType() == ReactionType::POLYMERIZATIONPLUSEND) {
                
                //If reaching a threshold time for manual treadmilling rate changer
                if(tau() > SysParams::DRParams.manualCharStartTime){
                    //all bare rate will be change by a threshold ratio
                    newRate = _polyChanger[_type]->changeRate(r->getBareRate() * SysParams::DRParams.manualPlusPolyRate, force);
                }
                else{
                    newRate = _polyChanger[_type]->changeRate(r->getBareRate(), force);
                }
                
                r->setRateScaled(newRate);
                r->updatePropensity();

            }
            
            //change all plus end depolymerization rates, not force dependent
            //If reaching a threshold time for manual treadmilling rate changer
            if(tau() > SysParams::DRParams.manualCharStartTime){
                if(r->getReactionType() == ReactionType::DEPOLYMERIZATIONPLUSEND) {
                    r->setRateScaled(r->getBareRate() * SysParams::DRParams.manualPlusDepolyRate);
                    r->updatePropensity();
                }
            }
        }
    
    }
    
    //load force from back (affects minus end polymerization)
    if(_minusEnd) {
        
        //get force of front bead
        force = _b1->getLoadForcesM();
        
        //change all plus end polymerization rates
        for(auto &r : _cCylinder->getInternalReactions()) {
            float newRate;
            if(r->getReactionType() == ReactionType::POLYMERIZATIONMINUSEND) {
                
                //If reaching a threshold time for manual treadmilling rate changer
                if(tau() > SysParams::DRParams.manualCharStartTime){
                    //all bare rate will be change by a threshold ratio
                    newRate = _polyChanger[_type]->changeRate(r->getBareRate() * SysParams::DRParams.manualMinusPolyRate, force);
                }
                else{
                    newRate = _polyChanger[_type]->changeRate(r->getBareRate(), force);
                }
                
                r->setRateScaled(newRate);
                r->updatePropensity();
                
            }
            
            //change all minus end depolymerization rates, not force dependent
            //If reaching a threshold time for manual treadmilling rate changer
            if(tau() > SysParams::DRParams.manualCharStartTime){

                if(r->getReactionType() == ReactionType::DEPOLYMERIZATIONMINUSEND) {
                    r->setRateScaled(r->getBareRate() * SysParams::DRParams.manualMinusDepolyRate);
                    r->updatePropensity();
                }
            }
        }

    }
}

bool Cylinder::isFullLength() {
    
#ifdef MECHANICS
    return areEqual(_mCylinder->getEqLength(), SysParams::Geometry().cylinderSize[_type]);
#else
    return true;
#endif
}

void Cylinder::printSelf() {
    
    cout << endl;
    
    cout << "Cylinder: ptr = " << this << endl;
    cout << "Cylinder ID = " << getId() << endl;
    cout << "Parent ptr = " << getParent() << endl;
    cout << "Coordinates = " << coordinate[0] << ", " << coordinate[1] << ", " << coordinate[2] << endl;
    
    if(_plusEnd) cout << "Is a plus end." << endl;
    if(_minusEnd) cout << "Is a minus end." << endl;
    
    if(_branchingCylinder != nullptr) cout << "Has a branching cylinder." << endl;
    
    cout << "Position = " << _position << endl;
    
    cout<< "Length "<<_mCylinder->getLength()<<endl;
    cout<< "Eq Length "<<_mCylinder->getEqLength()<<endl;
    cout<< "Eq Theta "<<_mCylinder->getEqTheta()<<endl;
    cout<<" Stretching constant "<<_mCylinder->getStretchingConst()<<endl;
    cout<<" Bending constant "<<_mCylinder->getBendingConst()<<endl;
    
    cout << endl;
    
#ifdef CHEMISTRY
    cout << "Chemical composition of cylinder:" << endl;
    _cCylinder->printCCylinder();
#endif
    
    cout << endl;
    
    cout << "Bead information..." << endl;
    
    _b1->printSelf();
    _b2->printSelf();
    
    cout << endl;
}

bool Cylinder::within(Cylinder* other, double dist) {
    
    //check midpoints
    if(twoPointDistancesquared(coordinate, other->coordinate) <= (dist * dist))
        return true;
    
    //briefly check endpoints of other
    if(twoPointDistancesquared(coordinate, other->_b1->vcoordinate()) <= (dist * dist) ||
       twoPointDistancesquared(coordinate, other->_b2->vcoordinate()) <= (dist * dist))
        return true;
    
    return false;
}

vector<FilamentRateChanger*> Cylinder::_polyChanger;
ChemManager* Cylinder::_chemManager = 0;
