
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

#ifndef MEDYAN_Output_h
#define MEDYAN_Output_h

#include <fstream>

#include "common.h"

#include "Parser.h"

///FORWARD DECLARATIONS
class CompartmentGrid;
class SubSystem;

/// To print a specified output into a file
/*!
 *  An output object, initialized by the Controller, can print a number of specific
 *  output formats, including current snapshot, forces, tensions, and birth times. 
 *  Upon destruction, the output file is closed.
 */

class Output {
protected:
    ofstream _outputFile; ///< The output file being used
    
    SubSystem* _subSystem = nullptr;
    
public:
    /// Constructor, which opens the output file
    Output(string outputFileName, SubSystem* s) {
        _outputFile.open(outputFileName);
        if(!_outputFile.is_open()) {
            cout << "There was an error opening file " << outputFileName
            << " for output. Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        cout << "Opening file " << outputFileName << endl;
        
        _subSystem = s;
    }
    /// Destructor, which closes the output file
    ~Output() {_outputFile.close();}
    
    /// To be implemented in sub classes
    virtual void print(int snapshot) = 0;
};

/// Print basic information about all Filament, Linker,
/// MotorGhost, and CaMKIIingPoint
class BasicSnapshot : public Output {

public:
    BasicSnapshot(string outputFileName, SubSystem* s) : Output(outputFileName, s) {}
    ~BasicSnapshot() {}
    
    virtual void print(int snapshot);
};

/// Print birth times of beads for each Filament, Linker,
/// MotorGhost, and CaMKIIingPoint
class BirthTimes : public Output {
    
public:
    BirthTimes(string outputFileName, SubSystem* s) : Output(outputFileName, s) {}
    ~BirthTimes() {}
    
    virtual void print(int snapshot);
};

/// Print forces on beads for each Filament
class Forces : public Output {
    
public:
    Forces(string outputFileName, SubSystem* s) : Output(outputFileName, s) {}
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
    Tensions(string outputFileName, SubSystem* s) : Output(outputFileName, s) {}
    ~Tensions() {}
    
    virtual void print(int snapshot);
};

/// Print wall tension for each pinned filament:
///                 k * l * nhat
/// where k is the stretching force constant of the pin, l is the current
/// vector distance away from the pin position for the pinned bead.
/// @note - nhat is a vector pointing from the direction of the boundary normal.
class WallTensions : public Output {
    
public:
    WallTensions(string outputFileName, SubSystem* s) :
                Output(outputFileName, s) {}
    ~WallTensions() {}
    
    
    virtual void print(int snapshot);
};

    
/// Print type of each species
class Types : public Output {
    
public:
    Types(string outputFileName, SubSystem* s) : Output(outputFileName, s) {}
    ~Types() {}
        
    virtual void print(int snapshot);
};


/// Print all chemical species in the system, including diffusing
/// and bulk species, filament, motors, linkers and camkiiers.
class Chemistry : public Output {

ChemistryData _chemData; ///< chemistry data of this system
CompartmentGrid* _grid; ///< compartment grid of the system
    
public:
    Chemistry(string outputFileName, SubSystem* s,
              ChemistryData chemData, CompartmentGrid* grid)
    
        : Output(outputFileName, s),
         _chemData(chemData), _grid(grid) {}
    
    ~Chemistry() {}

    virtual void print(int snapshot);
};


/// Print MotorGhost binding lifetimes
class MotorLifetimes : public Output {
    
public:
    MotorLifetimes(string outputFileName, SubSystem* s) : Output(outputFileName, s) {}
    ~MotorLifetimes() {}
    
    virtual void print(int snapshot);
};

/// Print MotorGhost walk lengths
class MotorWalkLengths : public Output {
    
public:
    MotorWalkLengths(string outputFileName, SubSystem* s) : Output(outputFileName, s) {}
    ~MotorWalkLengths() {}
    
    virtual void print(int snapshot);
};

/// Print Linker binding lifetimes
class LinkerLifetimes : public Output {
    
public:
    LinkerLifetimes(string outputFileName, SubSystem* s) : Output(outputFileName, s) {}
    ~LinkerLifetimes() {}
    
    virtual void print(int snapshot);
};

/// Print Filament turnover times
class FilamentTurnoverTimes : public Output {
    
public:
    FilamentTurnoverTimes(string outputFileName, SubSystem* s) : Output(outputFileName, s) {}
    ~FilamentTurnoverTimes() {}
    
    virtual void print(int snapshot);
};

#endif