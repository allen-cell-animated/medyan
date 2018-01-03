
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

#    include "Component.h"
#    include "GController.h"
#    include "SubSystem.h"

namespace test_public {

    class ComponentDummy: public Component {
    public:
        short type;

        ComponentDummy(short newType): type(newType) {}

        virtual void printSelf() override {
            cout << "This is a dummy component object with type "
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

        SysParams::GParams.cylinderNumMon.resize(1, 3);
    }

}

#  endif // Include guard
#endif // TESTING
