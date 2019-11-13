
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
#include "CaMKIIingPoint.h"
#include "Boundary.h"

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
		//cerr<<"========== CaMKIIBindingManager CallBack" <<endl; //Carlos verbose prints
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
		if(delta != 1 && delta != -1)
			cerr <<"========== CaMKIIBundlingManager CallBack ---- delta = " << delta <<endl; //Carlos verbose prints

		//update this cylinder
		Compartment* c = _cylinder->getCompartment();

		int count = Cylinder::getCylinders().size();

		for(auto &manager : c->getFilamentBindingManagers()) {

			if(dynamic_cast<CaMKIIBundlingManager*>(manager.get())) {

				CCylinder* cc = _cylinder->getCCylinder();

				//update binding sites
				if(delta == +1) {
//					int bef = manager->numBindingSites();
					manager->addPossibleBindings(cc, _bindingSite);
//					int aft = manager->numBindingSites();

//					cerr << "========MILLAD: UpdateCaMKIIerBundlingCallback ADDING  " << cc << "  " <<
//					(dynamic_cast<CaMKIICylinder*>(_cylinder) == NULL ? "normal" : "camkii_cylinder") << " count: " << count << "  before: " << bef << "  after: " << aft;

//					cerr << "---------------============----------" << endl;
//					for(auto ccc : Cylinder::getCylinders()) {
//						if(dynamic_cast<CaMKIICylinder*>(ccc)) {
//							cerr << "========MILLAD: UpdateCaMKIIerBundlingCallback List CCylinder: " << ccc->getCCylinder() << endl;
//						}
//					}
//					cerr << "---------------============----------" << endl;
				}

				else /* -1 */ {
//					int bef = manager->numBindingSites();
					manager->removePossibleBindings(cc, _bindingSite);
//					int aft = manager->numBindingSites();

//					cerr << "========MILLAD: UpdateCaMKIIerBundlingCallback REMOVING  " << cc << "  " <<
//						 (dynamic_cast<CaMKIICylinder*>(_cylinder) == NULL ? "normal" : "camkii_cylinder") << " count: " << count << "  before: " << bef << "  after: " << aft;


//					cerr << "---------------============----------" << endl;
//					for(auto ccc : Cylinder::getCylinders()) {
//						if(dynamic_cast<CaMKIICylinder*>(ccc)) {
//							cerr << "========MILLAD: UpdateCaMKIIerBundlingCallback List CCylinder: " << ccc->getCCylinder() << endl;
//						}
//					}
//					cerr << "---------------============----------" << endl;
				}
			}
		}

//		count = Cylinder::getCylinders().size();
//
//		cerr << "  UPDATED count: " << count << endl;

	}
};

struct UpdateCaMKIIerDummyCylinderCallback {

	Cylinder* _cylinder; ///< cylinder to update

	short _bindingSite;  ///< binding site to update

	//Constructor, sets members
	UpdateCaMKIIerDummyCylinderCallback(Cylinder* cylinder, short bindingSite)

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

				if(delta == -1 && n == 0) {
					manager->removePossibleBindings(cc, _bindingSite);
//					cerr << "====== MILLAD: Removing CaMKIICylinder from PossibleCylinders " << _cylinder << "  " << _cylinder->getCCylinder()
//						 << endl;
				}
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
//		cerr<<"========== CaMKII Unbinding CallBack "; //Carlos verbose prints
//		cerr<< "ID: " << _camkiiingPoint->getID() << " Coord:" << _camkiiingPoint->getCoordinationNumber() << endl; //Carlos verbose prints

		// Removing the internal reaction of the CaMKII cylinder
		_camkiiingPoint->getCaMKIICylinder()->getCCylinder()->removeInternalReaction(r);

		//remove the camkiiing point
		_camkiiingPoint->getCCaMKIIingPoint()->removeBond(_camkiiingPoint->getBond(0));
//		cerr << "====== MILLAD: Removing CaMKIICylinder " << _camkiiingPoint->getCaMKIICylinder()  << "   " << _camkiiingPoint->getCaMKIICylinder()->getCCylinder() << endl;
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
//		cerr<<"========== CaMKII Unbundling CallBack "; //Carlos verbose prints
//		cerr<< "ID: " << _camkiiingPoint->getID() << " Coord:" << _camkiiingPoint->getCoordinationNumber()<< "-->"; //Carlos verbose prints
		auto bond = _camkiiingPoint->removeRandomBond();
		// Adding the cylinder assigned to the removed bond?
//        _bManager->addPossibleBindings(get<0>(bond)->getCCylinder(), get<1>(bond));
		if(_camkiiingPoint->getCoordinationNumber() == 1L) {

			// passivate unbundling reaction
			r->passivateReaction();

			// Getting a handler to the current off-reaction
			ReactionBase *offRxnBinding = _camkiiingPoint->getCCaMKIIingPoint()->getOffRxnBinding();

			// Creating an off-reaction if the off-reaction does not exist
			// due to CaMKII moving between compartments
			if(offRxnBinding == nullptr) {
				_camkiiingPoint->getCCaMKIIingPoint()->createOffReactionBinding(nullptr, _ps);
				ReactionBase *offRxn = _camkiiingPoint->getCCaMKIIingPoint()->getOffRxnBinding();
				_camkiiingPoint->getCCaMKIIingPoint()->setOffReaction(offRxn);
				offRxnBinding = offRxn;
			}

			// Reactivating the unbinding reaction
			offRxnBinding->activateReaction();

			// Removing the current off-reaction from the internal list in CaMKII cylinder
//			get<0>(_camkiiingPoint->getBonds()[0])->getCCylinder()->removeInternalReaction(r);
			_camkiiingPoint->getCaMKIICylinder()->getCCylinder()->removeInternalReaction(r);
//            get<0>(_camkiiingPoint->getBonds()[0])->getCCylinder()->addInternalReaction(offRxnBinding);
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

	//short _plusEnd;        ///< Plus end marker of new cylinder

	float _onRate;         ///< Rate of the binding reaction
	float _offRate;        ///< Rate of the unbinding reaction

	CaMKIIBindingCallback(CaMKIIBindingManager* bManager, float onRate, float offRate, SubSystem* ps)
			: _ps(ps), _bManager(bManager), _onRate(onRate), _offRate(offRate) {}

	void operator() (ReactionBase *r) {
		CaMKIIingPoint *b;
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
//            << " Using default values for Arp2/3 complex - l=10.0nm, theta=70.7deg"
//            << endl;
			l = 10.0;
			t = M_PI_2;
		}
#else
		cout << "CaMKIIing initialization cannot occur unless mechanics is enabled. Using"
			<< " default values for Arp2/3 complex - l=10.0nm, theta=70.7deg"
			<< endl;
			double l = 10.0;
			double t = 1.22;
#endif
		double s = SysParams::Geometry().monomerSize[filType];

		// initial position of CaMKIIingPoint
		auto cp = camkiiProjection(n, p, l);

		if (SysParams::RUNSTATE == true) {

			b = _ps->addTrackable<CaMKIIingPoint>(c1, camkiiType, pos, cp);
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
				b = _ps->addTrackable<CaMKIIingPoint>(c1, camkiiType, pos, cp);

				x = c->getCMonomer(0);
				frate = 0.0;
			} else {
				cout << "CaMKIIer Error. Cannot find binding Site in the list. Cannot complete restart. Exiting."
					 << endl;
				exit(EXIT_FAILURE);
			}
		}

		//create off reaction
		auto cCaMKIIer = b->getCCaMKIIingPoint();

		cCaMKIIer->setRates(_onRate, frate);
		cCaMKIIer->createOffReactionCaMKII(r, _ps, nullptr);
		cCaMKIIer->getOffReaction()->setBareRate(SysParams::CaMKIIUnbindingBareRate[camkiiType]);
	}
};

/// Callback to create a CaMKIIingPoint on a Filament
struct CaMKIIBundlingCallback {

	SubSystem* _ps;        ///< ptr to subsystem

	CaMKIIBundlingManager* _bManager; ///< CaMKIIing manager for this compartment

	//short _plusEnd;        ///< Plus end marker of new cylinder

	float _onRate;         ///< Rate of the binding reaction
	float _offRate;        ///< Rate of the unbinding reaction

	CaMKIIBundlingCallback(CaMKIIBundlingManager* bManager, float onRate, float offRate, SubSystem* ps)
	//TODO add minSearchDist maxSearchDist maxcoordNumber as local data members
			: _ps(ps), _bManager(bManager), _onRate(onRate), _offRate(offRate) {}

	void operator() (ReactionBase *r) {
		CaMKIIingPoint* cp;
		// TODO CAMKII either find cp from NN or to get the one the reaction was on it
//        float frate;
		short camkiiType = _bManager->getBoundInt();

		//choose a random binding site from manager
		//TODO Find a CAMKII with coordination number less than max number
		auto site = _bManager->chooseBindingSites();

		Cylinder*       c1  = get<0>(get<0>(site))->getCylinder();
		short           pos1 = get<1>(get<0>(site));

		Cylinder*       c2  = get<0>(get<1>(site))->getCylinder();
		short           pos2 = get<1>(get<1>(site));

		bool b1 = (dynamic_cast<CaMKIICylinder*>(c1)!= nullptr);
		bool b2 = (dynamic_cast<CaMKIICylinder*>(c2)!= nullptr);

		// If none of them or both of them are CaMKIIPoint, do nothing!
		if((!b1 && !b2) || (b1 && b2))
			return;

		if(b1) {
			cp = dynamic_cast<CaMKIICylinder *>(c1)->getCaMKIIPointParent();

			// cp should be CaMKIIPoint
			cp->addBond(c2, pos2);
		} else {
			cp = dynamic_cast<CaMKIICylinder *>(c2)->getCaMKIIPointParent();

			// cp should be CaMKIIPoint
			cp->addBond(c1, pos1);
		}

//TODO CJY make sure isn't needed before cleaning
#if 0

		// CamKIIingPoint cylinder pair filament with

        //get info from site
        Cylinder* c1 = get<0>(site)->getCylinder();
        short filType = c1->getType();

        double pos = double(get<1>(site)) / SysParams::Geometry().cylinderNumMon[filType];
        if(SysParams::RUNSTATE==true){
        //Get a position and direction of a new filament
        auto x1 = c1->getFirstBead()->coordinate;
        auto x2 = c1->getSecondBead()->coordinate;

        //get original direction of cylinder
        auto p= midPointCoordinate(x1, x2, pos);
        vector<double> n = twoPointDirection(x1, x2);

        //get camkii projection
#ifdef MECHANICS
        //use mechanical parameters
        double l, t;
        if(SysParams::Mechanics().CaMKIIStretchingL.size() != 0) {
            l = SysParams::Mechanics().CaMKIIStretchingL[camkiiType];
            t = SysParams::Mechanics().CaMKIIBendingTheta[camkiiType];
        }
        else {
            // TODO
            cout << "CaMKIIing  bundling initialization cannot occur unless mechanical parameters are specified."
            << " Using default values for Arp2/3 complex - l=10.0nm, theta=70.7deg"
            << endl;
            l = 10.0;
            t = 1.22;
        }
#else
        cout << "CaMKIIing initialization cannot occur unless mechanics is enabled. Using"
        << " default values for Arp2/3 complex - l=10.0nm, theta=70.7deg"
        << endl;
        double l = 10.0;
        double t = 1.22;
#endif
        double s = SysParams::Geometry().monomerSize[filType];

        auto camkiiPosDir = camkiiProjection(n, p, l, s, t);
        auto bd = get<0>(camkiiPosDir); auto bp = get<1>(camkiiPosDir);

        //TODO Search for new filament
//        //create a new filament
//        Filament* f = _ps->addTrackable<Filament>(_ps, filType, bp, bd, true, true);
//
//        //mark first cylinder

//        auto bindingSite = _bManager->chooseBindingSite();
//        auto ccx = get<0>(bindingSite);
//        auto cylin = ccx->getCylinder();
//        short bs = get<1>(ccx);

//        Cylinder* c = f->getCylinderVector().front();
//        //c->getCCylinder()->getCMonomer(0)->speciesPlusEnd(_plusEnd)->up();
//
//        //create new camkii
//            CMonomer* x=c->getCCylinder()->getCMonomer(0);
//            for(auto p = 0; p <SysParams::Geometry().cylinderNumMon[filType];p++){
//                auto xx =  c->getCCylinder()->getCMonomer(p)->speciesBound(SysParams::Chemistry().camkiierBoundIndex[filType]);
//                auto yy =c->getCCylinder()->getCMonomer(p)->speciesCaMKIIer(camkiiType);
//                auto zz =c->getCCylinder()->getCMonomer(p)->speciesFilament(0);
//                //std::cout<<c->getID()<<" "<<p<<" "<<xx->getN()<<" "<<yy->getN()<<" "<<zz->getN()<<endl;
//                            }
//            std::cout<<x->speciesFilament(0)->getN()<<" "<<x->speciesMinusEnd(0)->getN()<<endl;
//            vector<Cylinder*> cy{c1,c};
        //TODO Update new CAMKII to update coordination number
        //TODO Create the bonds between CAMKII and filament
        //b= _ps->addTrackable<CaMKIIingPoint>(c1, camkiiType, pos);

            for(auto p = 0; p <SysParams::Geometry().cylinderNumMon[filType];p++){
                auto xx =  c->getCCylinder()->getCMonomer(p)->speciesBound(SysParams::Chemistry().camkiierBoundIndex[filType]);
                auto yy =c->getCCylinder()->getCMonomer(p)->speciesCaMKIIer(camkiiType);
                auto zz =c->getCCylinder()->getCMonomer(p)->speciesFilament(0);
                //std::cout<<c->getID()<<" "<<p<<" "<<xx->getN()<<" "<<yy->getN()<<" "<<zz->getN()<<endl;
            }
            //std::cout<<x->speciesFilament(0)->getN()<<" "<<x->speciesMinusEnd(0)->getN()<<endl;
        frate=_offRate;
        }
        else
        {
            CCylinder* c; auto check = false;
        vector<tuple<tuple<CCylinder*, short>, tuple<CCylinder*, short>>> BrT=_bManager->getbtuple();
            for(auto T:BrT){
                CCylinder* cx=get<0>(get<0>(T));
                double p = double(get<1>(get<0>(T)))/ double(SysParams::Geometry().cylinderNumMon[filType]);
                if(cx->getCylinder()->getID()==c1->getID() && p==pos){
                    c=get<0>(get<1>(T));
                    check = true;
                    break;
                }}
            if(check){
                auto cyl = c->getCylinder();
                //std::cout<<twoPointDistance(cyl->getFirstBead()->coordinate,cyl->getSecondBead()->coordinate)<<" ";
            //std::cout<<c->getCylinder()->getID()<<endl;
                CMonomer* x=c->getCMonomer(0);
                for(auto p = 0; p <SysParams::Geometry().cylinderNumMon[filType];p++){
                    auto xx =  c->getCMonomer(p)->speciesBound(SysParams::Chemistry().camkiierBoundIndex[filType]);
                    auto yy =c->getCMonomer(p)->speciesCaMKIIer(camkiiType);
                    auto zz =c->getCMonomer(p)->speciesFilament(0);
                    //std::cout<<c->getCylinder()->getID()<<" "<<p<<" "<<xx->getN()<<" "<<yy->getN()<<" "<<zz->getN()<<endl;
                }
                //std::cout<<x->speciesFilament(0)->getN()<<" "<<x->speciesMinusEnd(0)->getN()<<endl;
                //vector<Cylinder*> cy{c1,c->getCylinder()};
                b= _ps->addTrackable<CaMKIIingPoint>(c1, camkiiType, pos);

                x=c->getCMonomer(0);
                for(auto p = 0; p <SysParams::Geometry().cylinderNumMon[filType];p++){
                    auto xx =  c->getCMonomer(p)->speciesBound(SysParams::Chemistry().camkiierBoundIndex[filType]);
                    auto yy =c->getCMonomer(p)->speciesCaMKIIer(camkiiType);
                    auto zz =c->getCMonomer(p)->speciesFilament(0);
                    //std::cout<<c->getCylinder()->getID()<<" "<<p<<" "<<xx->getN()<<" "<<yy->getN()<<" "<<zz->getN()<<endl;
                }
                //std::cout<<x->speciesFilament(0)->getN()<<" "<<x->speciesMinusEnd(0)->getN()<<endl;
            frate=0.0;
            }
            else
                cout<<"CaMKIIer Error. Cannot find binding Site in the list. Cannot complete restart. Exiting." <<endl;
                //exit(EXIT_FAILURE);
        }
#endif

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
