
#ifndef MEDYAN_CAMKII_ChemCallbacks_h
#define MEDYAN_CAMKII_ChemCallbacks_h

#include "common.h"
#include "utility.h"

#include "SubSystem.h"
#include "CompartmentGrid.h"
#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"
#include "Linker.h"
#include "MotorGhost.h"
#include "BranchingPoint.h"
#include "Structure/CaMKII/CaMKIIingPoint.h"
#include "Boundary.h"
//#include "ReactionBase.h" //added by jl135

#include "../BindingManager.h"
#include "CaMKIIChemCallbacks.h"

#include "GController.h"
#include "MathFunctions.h"
#include "SysParams.h"
#include "Rand.h"

using namespace mathfunc;


#ifdef RSPECIES_SIGNALING

struct UpdateCaMKIIerBindingCallback {

	Cylinder* _cylinder; ///< cylinder to update

	short _bindingSite;  ///< binding site to update

	//Constructor, sets members
	UpdateCaMKIIerBindingCallback(Cylinder* cylinder, short bindingSite)

			: _cylinder(cylinder), _bindingSite(bindingSite) {}

	//callback
	void operator() (RSpecies *r, int delta) {
		//update this cylinder
		Compartment* c = _cylinder->getCompartment();

		for(auto &manager : c->getFilamentBindingManagers()) {

			if(dynamic_cast<CaMKIIBindingManager*>(manager.get())) {

				CCylinder* cc = _cylinder->getCCylinder();

				//update binding sites
				if(delta == +1) manager->addPossibleBindings(cc, _bindingSite);

				else /* -1 */manager->removePossibleBindings(cc, _bindingSite);
			}
		}
	}
};

struct UpdateCaMKIIerBundlingCallback {

	Cylinder* _cylinder; ///< cylinder to update

	short _bindingSite;  ///< binding site to update

	//Constructor, sets members
	UpdateCaMKIIerBundlingCallback(Cylinder* cylinder, short bindingSite)

			: _cylinder(cylinder), _bindingSite(bindingSite) {}

	//callback
	void operator() (RSpecies *r, int delta) {
		//update this cylinder
		Compartment* c = _cylinder->getCompartment();

		int count = Cylinder::getCylinders().size();

		for(auto &manager : c->getFilamentBindingManagers()) {

			if(dynamic_cast<CaMKIIBundlingManager*>(manager.get())) {

				CCylinder* cc = _cylinder->getCCylinder();

				//update binding sites
				if(delta == +1) manager->addPossibleBindings(cc, _bindingSite);

				else /* -1 */ manager->removePossibleBindings(cc, _bindingSite);
			}
		}

	}
};

struct UpdateCaMKIIerCylinderCallback {

	Cylinder* _cylinder; ///< cylinder to update

	short _bindingSite;  ///< binding site to update

	//Constructor, sets members
	UpdateCaMKIIerCylinderCallback(Cylinder* cylinder, short bindingSite)

			: _cylinder(cylinder), _bindingSite(bindingSite) {}

	//callback
	void operator() (RSpecies *r, int delta) {
		if(!dynamic_cast<CaMKIICylinder*>(_cylinder))
			return;

		const int n = r->getN();

		//update this cylinder
		Compartment* c = _cylinder->getCompartment();

		for(auto &manager : c->getFilamentBindingManagers()) {

			if(dynamic_cast<CaMKIIBundlingManager*>(manager.get())) {

				CCylinder* cc = _cylinder->getCCylinder();

				if(delta == -1 && n == 0)
					manager->removePossibleBindings(cc, _bindingSite);
			}
		}
	}
};

#endif


#ifdef REACTION_SIGNALING

/// Callback to unbind a CaMKIIingPoint from a Filament
struct CaMKIIingPointUnbindingCallback {

	SubSystem* _ps;
	CaMKIIingPoint* _camkiiingPoint;

	CaMKIIingPointUnbindingCallback(CaMKIIingPoint* b, SubSystem* ps)
			: _ps(ps), _camkiiingPoint(b) {}

	void operator() (ReactionBase *r) {
#ifdef DEBUG
		cout << "========== CaMKII Unbinding CallBack ID: " << _camkiiingPoint->getID() <<" Coord:" << _camkiiingPoint->getCoordinationNumber(); //Carlos verbose prints, jl135

#endif

		// Removing the internal reaction of the CaMKII cylinder
		_camkiiingPoint->getCaMKIICylinder()->getCCylinder()->removeInternalReaction(r);

		//remove the camkiiing point
		_camkiiingPoint->getCCaMKIIingPoint()->removeBond(_camkiiingPoint->getBond(0));
		_ps->removeTrackable<CaMKIICylinder>(_camkiiingPoint->getCaMKIICylinder());
		_ps->removeTrackable<CaMKIIingPoint>(_camkiiingPoint);
		delete _camkiiingPoint;
	}
};

/// Callback to unbind a CaMKIIingPoint from a Filament
struct CaMKIIingPointUnbundlingCallback {

	SubSystem* _ps;
	CaMKIIingPoint* _camkiiingPoint;
	CaMKIIBundlingManager* _bManager;

	CaMKIIingPointUnbundlingCallback(CaMKIIingPoint* b, SubSystem* ps, CaMKIIBundlingManager* bManager)
			: _ps(ps), _camkiiingPoint(b), _bManager(bManager) {}

	void operator() (ReactionBase *r) {
#ifdef DEBUG
		cout << "========== CaMKII Unbundling CallBack "; //Carlos verbose prints
		cout << "ID: " << _camkiiingPoint->getID() << " Coord:" << _camkiiingPoint->getCoordinationNumber() <<" CaMKII Type: " << _bManager->getBoundInt() <<" Max CaMKII Coord Number: " <<_bManager->getMaxCoordinationNumber() << endl; //Carlos verbose prints, jl135, check this
#endif

		auto bond = _camkiiingPoint->removeRandomBond(_bManager);

		// Adding the cylinder assigned to the removed bond?
		if(_camkiiingPoint->getCoordinationNumber() == 1L) {

			// passivate unbundling reaction
			r->passivateReaction();

			// Getting a handler to the current off-reaction
			ReactionBase *offRxnBinding = _camkiiingPoint->getCCaMKIIingPoint()->getOffRxnBinding();

			// Creating an off-reaction if the off-reaction does not exist
			// due to CaMKII moving between compartments
			if(offRxnBinding == nullptr) {
				_camkiiingPoint->getCCaMKIIingPoint()->setOffRate(_camkiiingPoint->getCCaMKIIingPoint()->getOffRateBinding());
				_camkiiingPoint->getCCaMKIIingPoint()->createOffReactionBinding(nullptr, _ps);
				ReactionBase *offRxn = _camkiiingPoint->getCCaMKIIingPoint()->getOffRxnBinding();
				_camkiiingPoint->getCCaMKIIingPoint()->setOffReaction(offRxn);
				offRxnBinding = offRxn;
			}

			// Reactivating the unbinding reaction
			offRxnBinding->activateReaction();

			// Removing the current off-reaction from the internal list in CaMKII cylinder
			_camkiiingPoint->getCaMKIICylinder()->getCCylinder()->removeInternalReaction(r);
			_camkiiingPoint->getCaMKIICylinder()->getCCylinder()->addInternalReaction(offRxnBinding);
			_camkiiingPoint->getCCaMKIIingPoint()->setOffReaction(offRxnBinding);

			// Removing the unbundling reaction
			_camkiiingPoint->getCCaMKIIingPoint()->setOffRxnBundling(nullptr);

		} else if(_camkiiingPoint->getCoordinationNumber() == _bManager->getMaxCoordinationNumber()-1 ) {
			_bManager->addPossibleBindings(_camkiiingPoint->getCaMKIICylinder()->getCCylinder());
		}

		_camkiiingPoint->updateReactionRates();
	}
};


/// Callback to create a CaMKIIingPoint on a Filament
struct CaMKIIBindingCallback {

	SubSystem* _ps;        ///< ptr to subsystem

	CaMKIIBindingManager* _bManager; ///< CaMKIIing manager for this compartment

	float _onRate;         ///< Rate of the binding reaction
	float _offRate;        ///< Rate of the unbinding reaction

	CaMKIIBindingCallback(CaMKIIBindingManager* bManager, float onRate, float offRate, SubSystem* ps)
			: _ps(ps), _bManager(bManager), _onRate(onRate), _offRate(offRate) {}

	void operator() (ReactionBase *r) {
		CaMKIIingPoint *camkii;
		float frate;
		short camkiiType = _bManager->getBoundInt();

		//choose a random binding site from manager
		auto site = _bManager->chooseBindingSite();
		//get info from site
		Cylinder *c1 = get<0>(site)->getCylinder();
		short filType = c1->getType();

		double pos = double(get<1>(site)) / SysParams::Geometry().cylinderNumMon[filType];

		//Get a position and direction of a new filament
		auto x1 = c1->getFirstBead()->coordinate;
		auto x2 = c1->getSecondBead()->coordinate;

		//get original direction of cylinder
		auto p = midPointCoordinate(x1, x2, pos);
		vector<double> n = twoPointDirection(x1, x2);

		//get camkii projection
#ifdef MECHANICS
		//use mechanical parameters
		double l, t;
		if (SysParams::Mechanics().CaMKIIStretchingL.size() != 0) {
			l = SysParams::Mechanics().CaMKIIStretchingL[camkiiType];
//            t = SysParams::Mechanics().CaMKIIBendingTheta[camkiiType];
			t = M_PI_2;
		} else {
//            cerr << "CaMKIIing initialization cannot occur unless mechanical parameters are specified."
//            << " Using default values for CaMKII complex - l=10.0nm, theta=90deg"
//            << endl;
			l = 10.0;
			t = M_PI_2;
		}
#else
		cerr << "CaMKIIing initialization cannot occur unless mechanical parameters are specified."
            << " Using default values for CaMKII complex - l=10.0nm, theta=90deg"
            << endl;
			double l = 10.0;
			double t = M_PI_2;
#endif
		double s = SysParams::Geometry().monomerSize[filType];

		// initial position of CaMKIIingPoint
		auto cp = camkiiProjection(n, p, l);

		if (SysParams::RUNSTATE == true) {

			camkii = _ps->addTrackable<CaMKIIingPoint>(c1, camkiiType, pos, cp);

#ifdef DEBUG
			cout << "========== CaMKII Binding CallBack - ID: " << camkii->getID() << endl; //Carlos verbose prints
#endif

			frate = _offRate;
		} else {
			CCylinder *c;
			auto check = false;
			vector<tuple<tuple<CCylinder *, short>, tuple<CCylinder *, short>>> BrT = _bManager->getbtuple();
			for (auto T:BrT) {
				CCylinder *cx = get<0>(get<0>(T));
				double p = double(get<1>(get<0>(T))) / double(SysParams::Geometry().cylinderNumMon[filType]);
				if (cx->getCylinder()->getID() == c1->getID() && p == pos) {
					c = get<0>(get<1>(T));
					check = true;
					break;
				}
			}
			if (check) {
				CMonomer *x = c->getCMonomer(0);
				vector<Cylinder *> cy{c1, c->getCylinder()};

				camkii = _ps->addTrackable<CaMKIIingPoint>(c1, camkiiType, pos, cp);

#ifdef DEBUG
				cout << "========== CaMKII Binding CallBack - ID: " << camkii->getID() <<" CaMKII Type: " << camkiiType << endl; //Carlos verbose prints, jl135
#endif

				x = c->getCMonomer(0);
				frate = 0.0;
			} else {
				cout << "CaMKIIer Error. Cannot find binding Site in the list. Cannot complete restart. Exiting."
					 << endl;
				exit(EXIT_FAILURE);
			}
		}

		//create off reaction
		auto cCaMKIIer = camkii->getCCaMKIIingPoint();

		cCaMKIIer->setOffRateBinding(frate);
		cCaMKIIer->setRates(_onRate, frate);
		cCaMKIIer->createOffReactionCaMKII(r, _ps, nullptr);
		cCaMKIIer->getOffReaction()->setBareRate(SysParams::CaMKIIUnbindingBareRate[camkiiType]);

	}
};

/// Callback to create a CaMKIIingPoint on a Filament
struct CaMKIIBundlingCallback {

	SubSystem* _ps;        ///< ptr to subsystem

	CaMKIIBundlingManager* _bManager; ///< CaMKIIing manager for this compartment

	float _onRate;         ///< Rate of the binding reaction
	float _offRate;        ///< Rate of the unbinding reaction

	CaMKIIBundlingCallback(CaMKIIBundlingManager* bManager, float onRate, float offRate, SubSystem* ps)
			: _ps(ps), _bManager(bManager), _onRate(onRate), _offRate(offRate) {}

	void operator() (ReactionBase *r) {

		short camkiiType = _bManager->getBoundInt();

		//choose a random binding site from manager
		auto site = _bManager->chooseBindingSites();

		Cylinder*       c1  = get<0>(get<0>(site))->getCylinder();
		short           pos1 = get<1>(get<0>(site));

		Cylinder*       c2  = get<0>(get<1>(site))->getCylinder();
		short           pos2 = get<1>(get<1>(site));

		bool b1 = (dynamic_cast<CaMKIICylinder*>(c1)!= nullptr);
		assert(b1 && "==== CaMKII Major Bug: the first site chosen by the bundling manager is not a CaMKIIPoint.");

		CaMKIIingPoint* cp;
		cp = dynamic_cast<CaMKIICylinder *>(c1)->getCaMKIIPointParent();

		//cout << "========== CaMKII Bundling CallBack without Flag  - "; //Carlos verbose prints
		//cout << "ID: " << cp->getID() << " Coord:" << cp->getCoordinationNumber() <<" CaMKII Type: " << camkiiType <<" Max CaMKII Coord Number: " <<_bManager->getMaxCoordinationNumber() << endl; //Carlos verbose prints, jl135, check this

		//r->getRNode()->printSelf(); //jl135
#ifdef DEBUG
		cout << "========== CaMKII Bundling CallBack - "; //Carlos verbose prints
		cout << "ID: " << cp->getID() << " Coord:" << cp->getCoordinationNumber() <<" CaMKII Type: " << camkiiType <<" Max CaMKII Coord Number: " <<_bManager->getMaxCoordinationNumber() << endl; //Carlos verbose prints, jl135, check this
#endif

		// cp should be CaMKIIPoint
		cp->addBond(c2, pos2);
		_bManager->removeAllCylindersOfFilamentFromCaMKIIPossibleBindings(c2, cp);

		//create off reaction
		auto cCaMKIIer = cp->getCCaMKIIingPoint();
		cCaMKIIer->setRates(_onRate, _offRate);
		cCaMKIIer->createOffReactionCaMKII(r, _ps, _bManager);
		cCaMKIIer->getOffReaction()->setBareRate(SysParams::CaMKIIUnbundlingBareRate[camkiiType]);
		_bManager->updateCaMKIIPossibleBindings(cp->getCaMKIICylinder()->getCCylinder());
	}
};

#endif

#endif //MEDYAN_CAMKII_ChemCallbacks_h
