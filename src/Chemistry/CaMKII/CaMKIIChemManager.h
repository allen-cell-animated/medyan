
#ifndef MEDYAN_CAMKII_ChemManager_h
#define MEDYAN_CAMKII_ChemManager_h

#include "common.h"

#include "ReactionTemplate.h"
#include "Parser.h"

//FORWARD DECLARATIONS
class Compartment;
class CompartmentGrid;
class ChemSim;

class Cylinder;
class CCylinder;
class CMonomer;


class CaMKIIChemManager {
public:
	static void setupBindingSites(ChemistryData &_chemData, int filType);
	static void setupBindingSitesInitCylinders(int filType);

	static void initCMonomerCaMKII(ChemistryData &_chemData, CMonomer* m, short filamentType, Compartment* c,
			int &bIndex);

	static void genFilBindingReactionsCaMKII(SubSystem* _subSystem, ChemistryData &_chemData, int filType,
														 CompartmentGrid* grid,
														 Compartment *C,
														 int &managerIndex, int &camkiiIndex);

	static void initializeManagers(SubSystem* _subSystem, Compartment* C0, int filType);

	static void genSpeciesCaMKII(ChemistryData &_chemData, Compartment& protoCompartment);

};

#endif //MEDYAN_CAMKII_ChemManager_h
