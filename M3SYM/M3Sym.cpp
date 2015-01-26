
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
 
 \subsection step3 Step 3: Running M3SYM
 
 See documentation on input files for more information. M3SYM
 executable must be run with the following command line arguments:
 
 -s : System input file to be used. Must be an absolute path
 -i : Input directory to be used, where all files specified in the
      system input file must be located. Must be an absolute path.
 -o : Output directory to be used (must be created beforehand),
      where all output files will be placed. Must be an absolute path.
 
 Run -h for help.
 
 */

#include <getopt.h>

#include "common.h"

#include "Controller.h"
#include "SubSystem.h"

void printUsage() {
    cout
    << "Usage: M3SYM -s systemFile -i inputDirectory -o outputDirectory"
    << endl;
}

int main(int argc, char **argv)
{
    
    cout << "******************** M3SYM **********************" << endl;
    
    cout.precision(8);
    
    //create subsystem and controller to run it
    SubSystem* s = nullptr;
    Controller c(s);

    string inputFile, inputDirectory, outputDirectory;
    int option;
    
    //parse command line args
    while ((option = getopt(argc, argv, "s:i:o:h")) != -1) {
        switch (option) {
            case 's' : inputFile = optarg;
                break;
            case 'i' : inputDirectory = optarg;
                break;
            case 'o' : outputDirectory = optarg;
                break;
            case 'h' : printUsage();
                exit(EXIT_FAILURE);
            default: printUsage();
                exit(EXIT_FAILURE);
        }
    }
    //check for arguments
    if(inputFile == "") {
        cout << "Must specify a system input file. Exiting." << endl;
        printUsage();
        exit(EXIT_FAILURE);
    }
    if(inputDirectory == "") {
        cout << "Must specify an input directory. Exiting." << endl;
        printUsage();
        exit(EXIT_FAILURE);
    }
    if(outputDirectory == "") {
        cout << "Must specify an output directory. Exiting." << endl;
        printUsage();
        exit(EXIT_FAILURE);
    }
    
    //initialize and run system
    c.initialize(inputFile, inputDirectory, outputDirectory);
    c.run();

}

