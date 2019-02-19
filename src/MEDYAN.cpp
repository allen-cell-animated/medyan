
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

/*! \mainpage MEDYAN software package
 
 \section intro_sec Introduction
 
The cell cytoskeleton plays a key role in human biology and disease, contributing ubiquitously
 to such important processes as embryonic development, wound repair and cancer 
 metastasis. Papoian laboratory is interested in gaining deeper understanding of the
 physical chemistry behind these complex, far-from-equilibrium mechanochemical 
 processes. Their approach and model, named Mechanochemical Dynamics of Active Networks
 (**MEDYAN**), is based on combining stochastic reaction-diffusion treatment
 of cellular biochemical processes with polymer physics of cytoskeletal filament network 
 growth, while explicitly coupling chemistry and mechanics.
 
 Papoian laboratory has developed a third-generation software package based on
 the **MEDYAN** model, to simulate growth dynamics of actin based filamentous networks *in vitro* 
 and *in vivo*. Recent papers where **MEDYAN** or its second-generation predecessor, **StochTools**,
 were used can be found on the publication section of [the Papoian group's main web page]
 (http://papoian.chem.umd.edu/ ) or on [the MEDYAN web page] (http://www.medyan.org ).
 The **MEDYAN** package can also be extended to simulate the dynamics of any active
 matter network.
 
 \section install_sec Installation
 
 \subsection step1 Step 1: Prerequisites
 
 The following libraries need to be installed first:
 See Installation guide (docs/InstallGuide.pdf) for more details.
 
 - Boost 1.49
 - GSL ...
 
 \subsection step2 Step 2: Installation of MEDYAN itself
 
 Untar the MEDYAN source code into some directory, enter
 into the "MEDYAN" folder and execute "make" from the command line.
 
 See Installation guide (docs/InstallGuide.pdf) for more information
 on setting command line compilation macros, compiler compatibility, etc.
 
 \subsection step3 Step 3: Running MEDYAN
 
 See the Usage guide (docs/UsageGuide.pdf) for more information. The MEDYAN
 executable must be run with the following command line arguments:
 
 -s : System input file to be used. Must be an absolute path
 -i : Input directory to be used, where all files specified in the
      system input file must be located. Must be an absolute path.
 -o : Output directory to be used (must be created beforehand),
      where all output files will be placed. Must be an absolute path.
 
 Run -h for help.
 
 */

#include "common.h"

#include "analysis/io/read_snapshot.h"
#include "core/controller/Controller.h"
#include "core/io/command_line.h"
#include "core/globals.h"
#include "Rand.h"
#include "Structure/SubSystem.h"
#include "utility.h"
#include "util/io/log.h"
#include "visual/window.hpp"

using namespace medyan;

int main(int argc, char **argv) {

    cout << endl;
    cout << "*********************** MEDYAN ************************" << endl;
    cout << "   Simulation package for the Mechanochemical Dynamics " << endl;
    cout << "         of Active Networks, Third Generation.         " << endl;
    cout << "         PAPOIAN LAB 2015, ALL RIGHTS RESERVED         " << endl;
    cout << "*******************************************************" << endl;
    
    cout.precision(8);
    
    // Parse command line (Will abort the program if parsing fails)
    commandline::readFromCommandLine(argc, argv);

    /**************************************************************************
    Initializations
    **************************************************************************/

    //create subsystem and controller to run it
    SubSystem* s = nullptr;
    Controller c(s);

    // Initialize the logger
    ::medyan::logger::Logger::defaultLoggerInitialization();

    // Seed global random generator
    if(!Global::readGlobal().randomGenSeedFixed) {
        Global::global().randomGenSeed = rdtsc();
        LOG(DEBUG) << "Global RNG seed: " << Global::readGlobal().randomGenSeed;
    }
    Rand::eng.seed(Global::readGlobal().randomGenSeed);

    // Visual
    visual::init();

    /**************************************************************************
    Start program 
    **************************************************************************/
    switch(Global::readGlobal().mode) {
    case GlobalVar::RunMode::Simulation:
        //initialize and run system
        c.initialize(Global::readGlobal().systemInputFile,
                     Global::readGlobal().inputDirectory,
                     Global::readGlobal().outputDirectory);
        c.run();
        break;
    case GlobalVar::RunMode::Analysis:
        {
            string inputFilePath = Global::readGlobal().inputDirectory + "/snapshot.traj";
            string pdbFilePath = Global::readGlobal().outputDirectory + "/snapshot.pdb";
            string psfFilePath = Global::readGlobal().outputDirectory + "/snapshot.psf";
            analysis::SnapshotReader sr(inputFilePath, pdbFilePath, psfFilePath);
            sr.readAndConvertToVmd();
        }
        break;
    }

}

