
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

#ifdef TESTING

#  ifndef MEDYAN_TEST_PUBLIC_H
#  define MEDYAN_TEST_PUBLIC_H

/* Some public classes and functions are defined here for test use.
*/

#    include "common.h"

#    include "CController.h"
#    include "Composite.h"
#    include "GController.h"
#    include "SubSystem.h"

namespace test_public {

    class CompositeDummy: public Composite {
    public:
        short type;

        CompositeDummy(short newType): type(newType) {}

        virtual void printSelf() override {
            cout << "This is a dummy composite object with type "
                 << type << endl;
        }
    
        virtual int getType() override { return type; }

    };

    // Initialize global variables to get a playground for test.
    // There is no need to undo this function to clean the test.
    void quickSetupPlayground(SubSystem* s, double size=1e10, size_t nx=1, size_t ny=1, size_t nz=1) {
        SysParams::GParams.compartmentSizeX = size;
        SysParams::GParams.compartmentSizeY = size;
        SysParams::GParams.compartmentSizeZ = size;
        
        SysParams::GParams.NX = nx;
        SysParams::GParams.NY = ny;
        SysParams::GParams.NZ = nz;
        
        SysParams::GParams.nDim = 3;

        GController g(s); // Dummy variable to initialize the compartments
        g.initializeGrid();
    }

    void quickSetupChem(SubSystem* s, string chemAlgorithm="GILLESPIE") {
        CController c(s);
        ChemistryData cData; // Dummy data
        c.initialize(chemAlgorithm, cData);
    }

}

#  endif // Include guard
#endif // TESTING
