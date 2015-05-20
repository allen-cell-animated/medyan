
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
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

/// For initailizing chemical reactions based on a specific system
/*!
 *  ChemManager is a singleton used for initailizing all chemistry in the system. 
 *  Initialized by the CController, the ChemManager initializes the chemical components 
 *  of the CompartmentGrid as well as initializes all [CCylinders] (@ref CCylinder) 
 *  created. The ChemManager can also update chemical components of a CCylinder.
 */
class ChemManager {
    
public:
    ///Set the chemManager instance
    static void setInstance(ChemManagerImpl *cii);
    
    ///Initialize the CompartmentGrid, based on the given simulation
    static void initializeSystem();
    
    ///Initializer for a Cylinder, based on the given simulation
    static void initializeCCylinder(CCylinder* cc, Filament* f,
                                    bool extensionFront,
                                    bool extensionBack,
                                    bool initialization);
    
    /// Update the copy numbers of all species in the chemical network
    /// @note - this only sets the copy number if the simulation time
    ///         tau has passed the release time of the molecule. This
    ///         function is called at every set of chemical steps to check
    ///         if molecules should be released at the current time.
    static void updateCopyNumbers();
    
private:
    static ChemManagerImpl* _pimpl; ///< Store a pointer to a specific implementation
                                    ///< of the initializer; no ownership
    ChemManager() {};

};


#endif