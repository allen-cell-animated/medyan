
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

/*! \mainpage M3SYM software package
 
 \section intro_sec Introduction
 
 Cell motility plays a key role in human biology and disease, contributing ubiquitously
 to such important processes as embryonic development, wound repair and cancer 
 metastasis. Papoian laboratory is interested in gaining deeper understanding of the
 physical chemistry behind these complex, far-from-equilibrium mechano-chemical 
 processes. His approach and model, named Mechano-chemical Dynamics of Active Networks,
 3rd Generation (MEDYAN3), is based on combining stochastic reaction-diffusion treatment
 of cellular biochemical processes with polymer physics of cytoskeletal filament network 
 growth, while explicitly coupling chemistry and mechanics.
 
 Papoian laboratory has developed **M3SYM**, a software package based on the MEDYAN3
 model, to simulate growth dynamics of actin based filamentous networks *in vitro* and 
 *in vivo*. Recent papers where **M3SYM** or its predecessor, **StochTools**, were used 
 can be found on the publication section of [the Papoian group's main web page: ]
 (http://papoian.chem.umd.edu/ ). The **M3SYM** package can also be extended to simulate
 the dynamics of any active matter network.
 
 \section install_sec Installation
 
 \subsection step1 Step 1: Prerequisites
 
 The following software packages need to be installed first:
 
 - Boost 1.49
 - GSL ...
 
 \subsection step2 Step 2: Installation of M3SYM itself
 
 Untar the **M3SYM** source code into some directory, enter 
 into the "M3SYM" and execute "make" from the command line.
 
 */

#include "common.h"

#include "Controller.h"
#include "SubSystem.h"

int main(int argc, const char * argv[])
{

    SubSystem s;
    Controller c(&s);

    c.initialize("/Users/jameskomianos/Code/M3SYM/M3SYM/",
                 "/Users/jameskomianos/Code/M3SYM/M3SYM/");
    c.run();

}

