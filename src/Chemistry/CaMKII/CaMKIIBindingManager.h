
#ifndef MEDYAN_CAMKII_BindingManager_h
#define MEDYAN_CAMKII_BindingManager_h

#include <unordered_set>
#include <unordered_map>
#include <random>

#include "common.h"

#include "Filament.h"
#include "NeighborListImpl.h"
#include "ReactionBase.h"

#include "SysParams.h"
#include "Rand.h"


/*
 * CJYM: We believe these constants point to the positions
 * of binding sites (bound species).
 *
 * For example, C1_RXN_INDEX is being used for CaMKII binding reactions.
 * The order of the reactant species for CaMKII binding reactions is defined
 * in ChemManager.cpp.
 * The index for the CaMKII binding site is 1.
 *
 */
#define C1_RXN_INDEX 1
#define C2_RXN_INDEX 0


/// Manager for Filament and CaMKIIingPoint creation
class CaMKIIBindingManager : public FilamentBindingManager {

	friend class ChemManager;

private:
	///Nucleation zone type, to define where nucleation should occur
	NucleationZoneType _nucleationZone;

	///If using a nucleation zone, nucleating distance from the boundary
	double _nucleationDistance;

	///possible bindings at current state
	unordered_set<tuple<CCylinder*, short>> _possibleBindings;
	vector<tuple<tuple<CCylinder*, short>, tuple<CCylinder*, short>>> _camkiirestarttuple; //Used only during restart conditions.
public:
	CaMKIIBindingManager(ReactionBase* reaction,
						 Compartment* compartment,
						 short boundInt, string boundName,
						 short filamentType,
						 NucleationZoneType zone = NucleationZoneType::ALL,
						 double nucleationDistance = numeric_limits<double>::infinity());
	~CaMKIIBindingManager() {}

	//@{
	///add possible binding reactions that could occur
	virtual void addPossibleBindings(CCylinder* cc, short bindingSite);
	virtual void addPossibleBindings(CCylinder* cc);
	//@}

	//@{
	/// Remove all bindings including this cylinder
	virtual void removePossibleBindings(CCylinder* cc, short bindingSite);
	virtual void removePossibleBindings(CCylinder* cc);
	//@}

	///update all possible binding reactions that could occur
	virtual void updateAllPossibleBindings();

	virtual int numBindingSites() {

		return _possibleBindings.size();
	}

	/// Choose a random binding site based on current state
	tuple<CCylinder*, short> chooseBindingSite() {

		assert((_possibleBindings.size() != 0)
			   && "Major bug: CaMKII Binding manager should not have zero binding \
                  sites when called to choose a binding site.");

		int randomIndex = Rand::randInteger(0, _possibleBindings.size() - 1);
		auto it = _possibleBindings.begin();

		advance(it, randomIndex);

		return *it;
	}

	virtual bool isConsistent();
	/// ARAVIND ADDED FEB 17 2016. append possible bindings.
	virtual void appendpossibleBindings(tuple<CCylinder*, short> t1, tuple<CCylinder*, short> t2){
		double oldN=numBindingSites();
		_possibleBindings.insert(t1);
		_camkiirestarttuple.push_back(make_tuple(t1,t2));
//        _camkiiCylinder=(get<0>(t2));
		double newN=numBindingSites();
		updateBindingReaction(oldN,newN);}
	virtual void clearpossibleBindings() {
		double oldN=numBindingSites();
		_possibleBindings.clear();
		updateBindingReaction(oldN,0);
	}
	vector<tuple<tuple<CCylinder*, short>, tuple<CCylinder*, short>>> getbtuple() {
		return _camkiirestarttuple;
	}
};

class CaMKIIBundlingManager : public FilamentBindingManager {

	friend class ChemManager;
	friend class CaMKIIChemManager;

private:
	float _rMin; ///< Minimum reaction range
	float _rMax; ///< Maximum reaction range
	size_t _maxCoordination;

	///possible bindings at current state
	unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>> _possibleBindings;
	vector<tuple<tuple<CCylinder*, short>, tuple<CCylinder*, short>>> _camkiirestarttuple; //Used only during restart conditions.

	//static neighbor list
	static vector<CylinderCylinderNL*> _neighborLists;
//    static vector<CaMKIIingPointCylinderNL*> _neighborLists;


	bool checkCoordinationNum(CCylinder *cc);

	// A function to find all the filaments that are connected to a CaMKII
	void getAllFilamentsOfCaMKII(CaMKIIingPoint *cp, unordered_set<Filament*> &s);

	void addPossibleBindingsNormal(CCylinder* normal_cc, short bindingSite, CaMKIIingPoint *inputCaMKII=nullptr);
	void addPossibleBindingsCaMKII(CCylinder* cc, short bindingSite);

	//Print binding sites, Aravind
	void printbindingsites();


public:
	CaMKIIBundlingManager(ReactionBase* reaction,
						  Compartment* compartment,
						  short boundInt,
						  string boundName,
						  short filamentType,
						  float rMax, float rMin, int maxCoordination);

	~CaMKIIBundlingManager() {}

	//@{
	///add possible binding reactions that could occur
	virtual void addPossibleBindings(CCylinder* cc, short bindingSite);
	virtual void addPossibleBindings(CCylinder* cc);
	//@}

	//@{
	/// Remove all bindings including this cylinder
	virtual void removePossibleBindings(CCylinder* cc, short bindingSite);
	virtual void removePossibleBindings(CCylinder* cc);
	//@}

	///update all possible binding reactions that could occur
	virtual void updateAllPossibleBindings();
	void updateCaMKIIPossibleBindings(CCylinder* cc);
	virtual int numBindingSites() {

		return _possibleBindings.size();
	}

	int checkLengthOfPossibleBindingInternal(unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>> &mm1,
											 unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>> &mm2);

	// This function checks whether possible bindings are the same before
	// and after updateAllPossibleBindings.
	void checkLengthOfPossibleBinding() {
		auto _poss_1st = _possibleBindings;
        auto len1=_possibleBindings.size();
        updateAllPossibleBindings();
        auto checkLength = checkLengthOfPossibleBindingInternal(_poss_1st, _possibleBindings);
        auto len2=_possibleBindings.size();

		cerr << "========== CaMKII checkLengthOfPossibleBinding - length: ";
        cerr << checkLength << "   len1: " << len1 << "   len2: " << len2 << endl;

		assert((len1 == len2)
			   && "Minor bug: CaMKIIing manager was not correctly updated when it reached this point (len1 == len2).");
		assert((len1 == checkLength)
			   && "Minor bug: CaMKIIing manager was not correctly updated when it reached this point (len1 == checkLength).");
	}

	/// Choose a random binding site based on current state
	tuple<tuple<CCylinder*, short>, tuple<CCylinder*, short>> chooseBindingSites();

	virtual bool isConsistent();
	/// ARAVIND ADDED FEB 17 2016. append possible bindings.
	virtual void appendpossibleBindings(tuple<CCylinder*, short> t1, tuple<CCylinder*, short> t2){
		// TODO fix restart for CAMKII
//		double oldN=numBindingSites();
		// _possibleBindings.emplace(t1,t2);
		// _camkiirestarttuple.push_back(make_tuple(t1,t2));
//        _camkiiCylinder=(get<0>(t2));
		//    double newN=numBindingSites();
		//   updateBindingReaction(oldN,newN);
	}
	virtual void clearpossibleBindings() {
		double oldN=numBindingSites();
		_possibleBindings.clear();
		updateBindingReaction(oldN,0);
	}

	vector<tuple<tuple<CCylinder*, short>, tuple<CCylinder*, short>>> getbtuple() {
		return _camkiirestarttuple;
	}

	//@{
	/// Getters for distances
	float getRMin() {return _rMin;}
	float getRMax() {return _rMax;}
	//@}

	//@{
	/// Getter for the max coordination number
	size_t getMaxCoordinationNumber() {return _maxCoordination;}
	//@}

	size_t getPossibleBindingSize() {
		return _possibleBindings.size();
	}



	/*
	 * Removing all cylinders of a filament from CaMKII PossibleBinding
	 */
	void removeAllCylindersOfFilamentFromCaMKIIPossibleBindings(Filament *f, CaMKIIingPoint *p);
	void removeAllCylindersOfFilamentFromCaMKIIPossibleBindings(Cylinder *c, CaMKIIingPoint *p);

	/*
	 * Adding all cylinders of a filament to CaMKII PossibleBinding
	 */
	void addAllCylindersOfFilamentToCaMKIIPossibleBindings(Filament *f, CaMKIIingPoint *p);
	void addAllCylindersOfFilamentToCaMKIIPossibleBindings(Cylinder *c, CaMKIIingPoint *p);
};

#endif //MEDYAN_CAMKII_BindingManager_h
