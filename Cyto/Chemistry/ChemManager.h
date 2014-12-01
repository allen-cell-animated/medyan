
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_ChemManager_h
#define M3SYM_ChemManager_h

#include "common.h"

#include "ChemManagerImpl.h"

//FORWARD DECLARATIONS
class Compartment;
class Filament;
class CCylinder;
struct ChemistryData;

/// ChemManager class is used for initailizing chemical reactions based on a specific system
/*!
 *  ChemManager is a singleton used for initailizing all chemistry in the system. Initialized by the [CController] (@ref CController),
 *  the ChemManager initializes the chemical components of the [CompartmentGrid] (@ref CompartmentGrid) as well as initializes all 
 *  [CCylinders] (@ref CCylinder) created. The ChemManager can also update chemical components of a [CCylinder] (@ref CCylinder).
 */
class ChemManager {
    
public:
    ///Set the chemManager instance
    static void setInstance(ChemManagerImpl *cii);
    
    ///Initialize the compartment grid, based on the given simulation
    static void initialize(ChemistryData& chem);
    
    ///Initializer, based on the given simulation
    static void initializeCCylinder(CCylinder* cc, Filament* f,bool extensionFront, bool extensionBack, bool creation);
    
    ///add/update cross cylinder reactions that are within range
    static void updateCCylinder(CCylinder* cc);
    
    
private:
    static ChemManagerImpl* _pimpl; ///< Store a pointer to a specific implementation of the initializer; no ownership
    ChemManager() {};

};


#endif