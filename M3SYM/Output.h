
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

#ifndef M3SYM_Output_h
#define M3SYM_Output_h

#include <fstream>

#include "common.h"

#include "Parser.h"

///FORWARD DECLARATIONS
class CompartmentGrid;

/// To print a specified output into a file
/*!
 *  An output object, initialized by the Controller, can print a number of specific
 *  output formats, including current snapshot, forces, tensions, and birth times. 
 *  Upon destruction, the output file is closed.
 */

class Output {
protected:
    ofstream _outputFile; ///< The output file being used
    
public:
    /// Constructor, which opens the output file
    Output(string outputFileName) {
        _outputFile.open(outputFileName);
        if(!_outputFile.is_open()) {
            cout << "There was an error opening file " << outputFileName
            << " for output. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        cout << "Opening file " << outputFileName << endl;
    }
    /// Destructor, which closes the output file
    ~Output() {_outputFile.close();}
    
    /// To be implemented in sub classes
    virtual void print(int step) = 0;
};

/// Print basic information about all Filament, Linker,
/// MotorGhost, and BranchingPoint
class BasicSnapshot : public Output {

public:
    BasicSnapshot(string outputFileName) : Output(outputFileName) {}
    ~BasicSnapshot() {}
    
    virtual void print(int step);
};

/// Print birth times of beads for each Filament, Linker,
/// MotorGhost, and BranchingPoint
class BirthTimes : public Output {
    
public:
    BirthTimes(string outputFileName) : Output(outputFileName) {}
    ~BirthTimes() {}
    
    virtual void print(int step);
};

/// Print forces on beads for each Filament
class Forces : public Output {
    
public:
    Forces(string outputFileName) : Output(outputFileName) {}
    ~Forces() {}
    
    virtual void print(int step);
};

/// Print tension for each Filament, Linker, and MotorGhost
/// @note This class prints the following:
///                 k * (l - l_0)
/// where k is the stretching force constant, l is the current
/// length, and l_0 is the equilibrium length.
class Tensions : public Output {
    
public:
    Tensions(string outputFileName) : Output(outputFileName) {}
    ~Tensions() {}
    
    virtual void print(int step);
};

/// Print all chemical species in the system, including diffusing
/// and bulk species, filament, motors, linkers and branchers.
class Chemistry : public Output {

ChemistryData _chemData; ///< chemistry data of this system
CompartmentGrid* _grid; ///< compartment grid of the system
    
public:
    Chemistry(string outputFileName,
              ChemistryData chemData,
              CompartmentGrid* grid)
    
        : Output(outputFileName),
         _chemData(chemData), _grid(grid) {}
    
    ~Chemistry() {}

    virtual void print(int step);
};

#endif
