
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

#include "Cylinder.h"

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
// beads but are free to be reallocated now.
void Cylinder::revectorize(cylinder* cylindervec, Cylinder** cylinderpointervec,
                        CCylinder** ccylindervec){
	#ifdef CROSSCHECK
	cout<<"revectorize"<<endl;
	#endif
    int i = 0;
    for(auto cyl:_cylinders.getElements()){
        //set _dcIndex
        cyl->_dcIndex = i;
        //copy attributes to a structure array
        cylindervec[i].filamentID = dynamic_cast<Filament*>(cyl->getParent())->getID();
        cylindervec[i].filamentposition = cyl->getPosition();
        cylindervec[i].bindices[0] = cyl->getFirstBead()->_dbIndex;
        cylindervec[i].bindices[1] = cyl->getSecondBead()->_dbIndex;
        cylindervec[i].cmpID = cyl->getCompartment()->getID();
        cylindervec[i].cindex = i;
        auto coord = cyl->coordinate;
        cylindervec[i].coord[0] = coord[0];
        cylindervec[i].coord[1] = coord[1];
        cylindervec[i].coord[2] = coord[2];
        cylindervec[i].type = cyl->getType();
        cylindervec[i].ID = cyl->getID();
        //other arrays needed
        ccylindervec[i] = cyl->getCCylinder();
        cylinderpointervec[i] = cyl;
        i++;
    }
    removedcindex.clear();
    Ncyl = _cylinders.getElements().size();
    maxcindex = Ncyl;
}

void Cylinder::appendrevectorize(cylinder* cylindervec, Cylinder** cylinderpointervec,
                           CCylinder** ccylindervec){
	cout<<"append revectorize"<<endl;
	int i = 0;
	maxcindex = 0;
	for(auto cyl:_cylinders.getElements()){

		i = cyl->_dcIndex;
		maxcindex = max(maxcindex,i);
		//copy attributes to a structure array
		cylindervec[i].filamentID = dynamic_cast<Filament*>(cyl->getParent())->getID();
		cylindervec[i].filamentposition = cyl->getPosition();
		cylindervec[i].bindices[0] = cyl->getFirstBead()->_dbIndex;
		cylindervec[i].bindices[1] = cyl->getSecondBead()->_dbIndex;
		cylindervec[i].cmpID = cyl->getCompartment()->getID();
		cylindervec[i].cindex = i;
		auto coord = cyl->coordinate;
		cylindervec[i].coord[0] = coord[0];
		cylindervec[i].coord[1] = coord[1];
		cylindervec[i].coord[2] = coord[2];
		cylindervec[i].type = cyl->getType();
		cylindervec[i].ID = cyl->getID();
		//other arrays needed
		ccylindervec[i] = cyl->getCCylinder();
		cylinderpointervec[i] = cyl;

	}
	maxcindex++;
	//Remove cindices in removedcindex vector that are greater-than-or-equal-to maxcindex
	for(auto reusablecidx=removedcindex.begin();reusablecidx != removedcindex.end();){
		if(*reusablecidx >= maxcindex)
			removedcindex.erase(reusablecidx);
		else
			reusablecidx++;
	}


	Ncyl = _cylinders.getElements().size();
}

void  Cylinder::copytoarrays() {
    long i =_dcIndex;
    cylinder* cylindervec = CUDAcommon::serlvars.cylindervec;
    Cylinder** cylinderpointervec = CUDAcommon::serlvars.cylinderpointervec;
    CCylinder** ccylindervec = CUDAcommon::serlvars.ccylindervec;
    //copy attributes to a structure array
    cylindervec[i].filamentID = dynamic_cast<Filament*>(this->getParent())->getID();
    cylindervec[i].filamentposition = _position;
    cylindervec[i].bindices[0] = _b1->_dbIndex;
    cylindervec[i].bindices[1] = _b2->_dbIndex;
    cylindervec[i].cmpID = _compartment->getID();
    cylindervec[i].cindex = i;
    cylindervec[i].type = _type;
    cylindervec[i].ID = _ID;
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
    Cylinder** cylinderpointervec = CUDAcommon::serlvars.cylinderpointervec;
    CCylinder** ccylindervec = CUDAcommon::serlvars.ccylindervec;
    resetcylinderstruct(cylindervec, _dcIndex);
    cylinderpointervec[_dcIndex] = NULL;
    ccylindervec[_dcIndex] = NULL;
}

void Cylinder::updateCoordinate() {
    coordinate = midPointCoordinate(_b1->coordinate, _b2->coordinate, 0.5);
    //update the coordiante in cylinder structure.
    cylinder* cylindervec = CUDAcommon::serlvars.cylindervec;
    cylindervec[_dcIndex].coord[0] = coordinate[0];
    cylindervec[_dcIndex].coord[1] = coordinate[1];
    cylindervec[_dcIndex].coord[2] = coordinate[2];
}


Cylinder::Cylinder(Composite* parent, Bead* b1, Bead* b2, short type, int position,
                   bool extensionFront, bool extensionBack, bool initialization)

    : Trackable(true, true, true, false),
      _b1(b1), _b2(b2), _type(type), _position(position), _ID(_cylinders.getID()) {

	//@{

    parent->addChild(unique_ptr<Component>(this));

    //revectorize if needed
    revectorizeifneeded();
    //set cindex
/*	_dcIndex = maxcindex;
	maxcindex++;*/
	//The following protocol is commented as it leads to seg faults.
	//Binding sites are stored based on cIndices and if a cIndex were to be reassigned
	// during chemistry, all the entries in the bindingmanager map go out of use. Hence,
	// it is necessary to not use this during chemistry.
    //set cindex based on maxbindex if there were no cylinders removed.
    if(removedcindex.size() == 0)
    {_dcIndex = maxcindex;
        maxcindex++;
    }
        // if cylinders were removed earlier, allot one of the available bead indices.
    else{
        _dcIndex = *removedcindex.begin();
        removedcindex.erase(removedcindex.begin());
    }
#ifdef CROSSCHECK
	cout<<"cindex "<<_dcIndex<< " alloted to ID "<<_ID<<endl;
#endif
    //@}

    //@{
    Ncyl = _cylinders.getElements().size();
    //check if you need to revectorize.
    cylinder* cylindervec = CUDAcommon::serlvars.cylindervec;
    Cylinder** cylinderpointervec = CUDAcommon::serlvars.cylinderpointervec;
    CCylinder** ccylindervec = CUDAcommon::serlvars.ccylindervec;
    //copy attributes to a structure array

    cylindervec[_dcIndex].filamentID = static_cast<Filament*>(this->getParent())->getID();
    cylindervec[_dcIndex].filamentposition = _position;
    cylindervec[_dcIndex].bindices[0] = _b1->_dbIndex;
    cylindervec[_dcIndex].bindices[1] = _b2->_dbIndex;

    //Set coordinate
    updateCoordinate();

    try {_compartment = GController::getCompartment(coordinate);}
    catch (exception& e) {
        cout << e.what() << endl;
        exit(EXIT_FAILURE);
    }

   //add to compartment
   _compartment->addCylinder(this);

    //@}
#ifdef MECHANICS
          //set eqLength according to cylinder size
              floatingpoint eqLength  = twoPointDistance(b1->coordinate, b2->coordinate);
          if(!SysParams::RUNSTATE) //RESTARTPHASE
          {
              int nummonomers = (int) round(eqLength/ SysParams::Geometry().monomerSize[type]);
              floatingpoint tpd = eqLength;
              
              if(nummonomers ==0){
                  eqLength = SysParams::Geometry().monomerSize[type];
              }
              else{
                  eqLength = (nummonomers) * SysParams::Geometry().monomerSize[type];
                  floatingpoint mindis = abs(tpd - eqLength);

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
    cylindervec[_dcIndex].cmpID = _compartment->getID();
    cylindervec[_dcIndex].cindex = _dcIndex;
    cylindervec[_dcIndex].type = _type;
    cylindervec[_dcIndex].ID = _ID;
    //other arrays needed
    ccylindervec[_dcIndex] = _cCylinder.get();
    cylinderpointervec[_dcIndex] = this;

    //init using chem manager
    _chemManager->initializeCCylinder(_cCylinder.get(), extensionFront,
                                      extensionBack, initialization);
#endif
    /*//copy further components to the array
    cylindervec[_dcIndex].cmpID = _compartment->getID();
    cylindervec[_dcIndex].cindex = _dcIndex;
    cylindervec[_dcIndex].type = _type;
    cylindervec[_dcIndex].ID = _ID;
    //other arrays needed
    ccylindervec[_dcIndex] = _cCylinder.get();
    cylinderpointervec[_dcIndex] = this;*/
}

Cylinder::~Cylinder() noexcept {
    //remove from compartment
    _compartment->removeCylinder(this);
    
}

/// Get filament type
int Cylinder::getType() {return _type;}

void Cylinder::updatePosition() {
	if(!setpositionupdatedstate) {

		//check if were still in same compartment, set new position
		updateCoordinate();

		Compartment *c;
		try { c = GController::getCompartment(coordinate); }
		catch (exception &e) {
			cout << e.what();

			printSelf();

			exit(EXIT_FAILURE);
		}

		if (c != _compartment) {
			mins = chrono::high_resolution_clock::now();

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
			//Removed March 8, 2019 Aravind. Unnecessary as all UpdatePosition calls are
			// immediately followed by UpdateNeighborLists call in Controller.cpp/.cu
/*        for(auto &manager : oldCompartment->getFilamentBindingManagers()) {
#ifdef NLORIGINAL
            manager->removePossibleBindings(oldCCylinder);
#endif
#ifdef NLSTENCILLIST
            manager->removePossibleBindingsstencil(oldCCylinder);
#endif
        }*/

			//clone and set new ccylinder
			CCylinder *clone = _cCylinder->clone(c);
			setCCylinder(clone);

			auto newCCylinder = _cCylinder.get();

        std::cout<<"moving cylinder with cindex "<<_dcIndex<<" and ID "<<_ID<<endl;
			//change both CCylinder and Compartment ID in the vector
			CUDAcommon::serlvars.cylindervec[_dcIndex].cmpID = _compartment->getID();
			CUDAcommon::serlvars.cylinderpointervec[_dcIndex] = this;
			CUDAcommon::serlvars.ccylindervec[_dcIndex] = _cCylinder.get();

			cout<<"Done "<<endl;
			//Add new ccylinder to binding managers
/*        for(auto &manager : newCompartment->getFilamentBindingManagers()){
#ifdef NLORIGINAL
            manager->addPossibleBindings(newCCylinder);
#endif
#ifdef NLSTENCILLIST
            //This directs call to Hybrid Binding Manager.
            manager->addPossibleBindingsstencil(newCCylinder);
#endif
        }*/
			mine = chrono::high_resolution_clock::now();
			chrono::duration<floatingpoint> compartment_update(mine - mins);
			CUDAcommon::tmin.timecylinderupdate += compartment_update.count();
			CUDAcommon::tmin.callscylinderupdate++;
		}
#endif

#ifdef MECHANICS
		//update length
		_mCylinder->setLength(twoPointDistance(_b1->coordinate,
		                                       _b2->coordinate));
#endif
	}
}

/// @note -  The function uses the bead load force to calculate this changed rate.
/// If there is no force on the beads the reaction rates are set to the bare.

void Cylinder::updateReactionRates() {
    
    floatingpoint force;
    
    //if no rate changer was defined, skip
    if(_polyChanger.empty()) return;
    
    //load force from front (affects plus end polymerization)
    if(_plusEnd) {
        
        //get force of front bead
        force = _b2->getLoadForcesP();
        
        //change all plus end polymerization rates
        for(auto &r : _cCylinder->getInternalReactions()) {
            
            if(r->getReactionType() == ReactionType::POLYMERIZATIONPLUSEND) {
            
                float newRate = _polyChanger[_type]->changeRate(r->getBareRate(), force);
                
                r->setRate(newRate);
                r->updatePropensity();
            }
        }
    }
    
    //load force from back (affects minus end polymerization)
    if(_minusEnd) {
        
        //get force of front bead
        force = _b1->getLoadForcesM();
        
        //change all plus end polymerization rates
        for(auto &r : _cCylinder->getInternalReactions()) {
            
            if(r->getReactionType() == ReactionType::POLYMERIZATIONMINUSEND) {
                
                float newRate =  _polyChanger[_type]->changeRate(r->getBareRate(), force);
                
                r->setRate(newRate);
                r->updatePropensity();
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
    cout << "Cylinder ID = " << _ID << endl;
    cout << "Parent ptr = " << getParent() << endl;
    cout << "Coordinates = " << coordinate[0] << ", " << coordinate[1] << ", " << coordinate[2] << endl;
    
    if(_plusEnd) cout << "Is a plus end." << endl;
    if(_minusEnd) cout << "Is a minus end." << endl;
    
    if(_branchingCylinder != nullptr) cout << "Has a branching cylinder." << endl;
    
    cout << "Position = " << _position << endl;
    
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

bool Cylinder::within(Cylinder* other, floatingpoint dist) {
    
    //check midpoints
    if(twoPointDistance(coordinate, other->coordinate) <= dist)
        return true;
    
    //briefly check endpoints of other
    if(twoPointDistance(coordinate, other->_b1->coordinate) <= dist ||
       twoPointDistance(coordinate, other->_b2->coordinate) <= dist)
        return true;
    
    return false;
}

vector<FilamentRateChanger*> Cylinder::_polyChanger;
ChemManager* Cylinder::_chemManager = 0;

Database<Cylinder*> Cylinder::_cylinders;
bool Cylinder::setpositionupdatedstate = false;
floatingpoint Cylinder::timecylinder1 = 0.0;
floatingpoint Cylinder::timecylinder2= 0.0;
floatingpoint Cylinder::timecylinderchem= 0.0;
floatingpoint Cylinder::timecylindermech= 0.0;
