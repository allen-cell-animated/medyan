
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
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
#include "DissipationTracker.h"


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
/// MotorGhost, and BranchingPoint
class BasicSnapshot : public Output {

public:
    BasicSnapshot(string outputFileName, SubSystem* s) : Output(outputFileName, s) {}
    ~BasicSnapshot() {}

    virtual void print(int snapshot);
};

/// Print birth times of beads for each Filament, Linker,
/// MotorGhost, and BranchingPoint
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
/// and bulk species, filament, motors, linkers and branchers.
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

/// Print total, chemdiss, mechdiss, chem, and mech
class Dissipation : public Output {

    ChemSim* _cs;

public:
    Dissipation(string outputFileName, SubSystem* s, ChemSim* cs)

    : Output(outputFileName, s), _cs(cs) {}

    ~Dissipation() {}

    virtual void print(int snapshot);
};

/// Print Filament plusend types
class PlusEnd : public Output {

public:
    PlusEnd(string outputFileName, SubSystem* s): Output(outputFileName, s) {}
    ~PlusEnd() {}

    virtual void print(int snapshot);
};


/// Print reactions for each filament
class ReactionOut : public Output {

public:
    ReactionOut(string outputFileName, SubSystem* s): Output(outputFileName, s) {}
    ~ReactionOut() {}
    
    virtual void print(int snapshot);
};


/// Print chem energy changes by HRCDID
class HRCD : public Output {

    ChemSim* _cs;

public:
    HRCD(string outputFileName, SubSystem* s, ChemSim* cs)

    : Output(outputFileName, s), _cs(cs) {}

    ~HRCD() {}

    virtual void print(int snapshot);
};


/// Print mech energy changes by HRMD
class HRMD : public Output {
    
    ChemSim* _cs;
    
public:
    HRMD(string outputFileName, SubSystem* s, ChemSim* cs)
    
    : Output(outputFileName, s), _cs(cs) {}
    
    ~HRMD() {}
    
    virtual void print(int snapshot);
};




// Print cm graph
class CMGraph : public Output {

public:
    CMGraph(string outputFileName, SubSystem* s)

    : Output(outputFileName, s) {}

    ~CMGraph() {}

    virtual void print(int snapshot);
};


// Print tm graph
class TMGraph : public Output {
    
public:
    TMGraph(string outputFileName, SubSystem* s)
    
    : Output(outputFileName, s) {}
    
    ~TMGraph() {}
    
    virtual void print(int snapshot);
};



// Print boundary repulsion force
class BRForces : public Output {

public:
    BRForces(string outputFileName, SubSystem* s) : Output(outputFileName, s) {}
    ~BRForces() {}

    virtual void print(int snapshot);
};


// Print concentration in each compartment
class Concentrations : public Output {

    SubSystem* _subSystem;///< SubSystem ptr
    ChemistryData _chemData; ///< chemistry data of this system

public:
    Concentrations(string outputFileName, SubSystem* s,
                   ChemistryData chemData)
    : Output(outputFileName, s), _subSystem(s), _chemData(chemData) {}
    ~Concentrations() {}

    virtual void print(int snapshot);
};

/// Print total, chemdiss, mechdiss, chem, and mech
class MotorWalkingEvents : public Output {
    
    ChemSim* _cs;
    
public:
    MotorWalkingEvents(string outputFileName, SubSystem* s, ChemSim* cs)
    
    : Output(outputFileName, s), _cs(cs) {}
    
    ~MotorWalkingEvents() {}
    
    virtual void print(int snapshot);
};

class LinkerUnbindingEvents : public Output {
    
    ChemSim* _cs;
    
public:
    LinkerUnbindingEvents(string outputFileName, SubSystem* s, ChemSim* cs)
    
    : Output(outputFileName, s), _cs(cs) {}
    
    ~LinkerUnbindingEvents() {}
    
    virtual void print(int snapshot);
};


class LinkerBindingEvents : public Output {
    
    ChemSim* _cs;
    
public:
    LinkerBindingEvents(string outputFileName, SubSystem* s, ChemSim* cs)
    
    : Output(outputFileName, s), _cs(cs) {}
    
    ~LinkerBindingEvents() {}
    
    virtual void print(int snapshot);
};



class HessianMatrix : public Output {
    
    ForceFieldManager* _ffm;
    
public:
    HessianMatrix(string outputFileName, SubSystem* s, ForceFieldManager* ffm)
    
    : Output(outputFileName, s), _ffm(ffm) {}
    
    ~HessianMatrix() {}
    
    virtual void print(int snapshot);
};


class HessianSpectra : public Output {
    
    ForceFieldManager* _ffm;
    
public:
    HessianSpectra(string outputFileName, SubSystem* s, ForceFieldManager* ffm)
    
    : Output(outputFileName, s), _ffm(ffm) {}
    
    ~HessianSpectra() {}
    
    virtual void print(int snapshot);
};





class Datadump : public Output {

    ChemSim* _cs;

public:
    Datadump(string outputFileName, SubSystem* s, ChemistryData chemData)

            : Output(outputFileName, s), _chemData(chemData),
            _outputFileName(outputFileName){}

    ~Datadump() {}

	ChemistryData _chemData; ///< chemistry data of this system
	string _outputFileName;

public:
    virtual void print(int snapshot);
};

#endif
