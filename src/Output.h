
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.0
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_Output_h
#define MEDYAN_Output_h

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
    virtual void print(int snapshot) = 0;
};

/// Print basic information about all Filament, Linker,
/// MotorGhost, and BranchingPoint
class BasicSnapshot : public Output {

public:
    BasicSnapshot(string outputFileName) : Output(outputFileName) {}
    ~BasicSnapshot() {}
    
    virtual void print(int snapshot);
};

/// Print birth times of beads for each Filament, Linker,
/// MotorGhost, and BranchingPoint
class BirthTimes : public Output {
    
public:
    BirthTimes(string outputFileName) : Output(outputFileName) {}
    ~BirthTimes() {}
    
    virtual void print(int snapshot);
};

/// Print forces on beads for each Filament
class Forces : public Output {
    
public:
    Forces(string outputFileName) : Output(outputFileName) {}
    ~Forces() {}
    
    virtual void print(int snapshot);
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
    
    virtual void print(int snapshot);
};

// Qin, 07/05/2016
class PlusEnd : public Output {
    
public:
    PlusEnd(string outputFileName): Output(outputFileName) {}
    ~PlusEnd() {}
    
    virtual void print(int snapshot);
};

/// Print all chemical species in the system, including diffusing
/// and bulk species, filament, motors, linkers and branchers.
class Chemistry : public Output {

ChemistryData _chemData; ///< chemistry data of this system
CompartmentGrid* _grid; ///< compartment grid of the system
    
public:
    Chemistry(string outputFileName, ChemistryData chemData,
                                     CompartmentGrid* grid)
    
        : Output(outputFileName),
         _chemData(chemData), _grid(grid) {}
    
    ~Chemistry() {}

    virtual void print(int snapshot);
};


/// Print MotorGhost binding lifetimes
class MotorLifetimes : public Output {
    
public:
    MotorLifetimes(string outputFileName) : Output(outputFileName) {}
    ~MotorLifetimes() {}
    
    virtual void print(int snapshot);
};

/// Print MotorGhost walk lengths
class MotorWalkLengths : public Output {
    
public:
    MotorWalkLengths(string outputFileName) : Output(outputFileName) {}
    ~MotorWalkLengths() {}
    
    virtual void print(int snapshot);
};

/// Print Linker binding lifetimes
class LinkerLifetimes : public Output {
    
public:
    LinkerLifetimes(string outputFileName) : Output(outputFileName) {}
    ~LinkerLifetimes() {}
    
    virtual void print(int snapshot);
};

/// Print Filament turnover times
class FilamentTurnoverTimes : public Output {
    
public:
    FilamentTurnoverTimes(string outputFileName) : Output(outputFileName) {}
    ~FilamentTurnoverTimes() {}
    
    virtual void print(int snapshot);
};

#endif
