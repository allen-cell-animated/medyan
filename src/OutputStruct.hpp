#ifndef MEDYAN_OutputStruct_Hpp
#define MEDYAN_OutputStruct_Hpp

#include <array>
#include <iostream>
#include <sstream>
#include <tuple>
#include <vector>

#include "MathFunctions.h"
#include "Structure/SurfaceMesh/Membrane.hpp"

///FORWARD DECLARATIONS
class BranchingPoint;
class Bubble;
class Filament;
class Linker;
class Membrane;
class MotorGhost;

/*!
 * The OutputStruct class serves as a connection between the real structure in
 * MEDYAN and the output structure. It may contain minimalistic data for
 * converting from/to system data and converting from/to output literals.
 * 
 * How much required information of each column can be obtained from each row
 * 
 *              | System Data   | Stored Data   | Output Data
 * -------------|---------------|---------------|---------------
 * System Data  | --            | Full          | Full
 * Stored Data  | Partial (#)   | --            | Full
 * Output Data  | Partial (#)   | Full          | --
 * 
 * (Except for snapshot serial which is not contained in the system)
 * 
 * (#): Not implemented
 */

class OutputStruct {
public:
    virtual ~OutputStruct() = default;

    /// get data from system data
    virtual void getFromSystem() { getFromSystemWithoutChildren(); }
    virtual void getFromSystemWithoutChildren() = 0; ///< Used in outputFromSystem

    /// output data from system data, with or without storing it
    /// This function SHOULD NOT call getFromSystem or outputFromStored inside.
    virtual void outputFromSystem(std::ostream& os) = 0;

    /// output stored data
    virtual void outputFromStored(std::ostream& os) { outputFromStoredWithoutChildren(os); }
    virtual void outputFromStoredWithoutChildren(std::ostream& os) = 0; ///< Used in outputFromSystem

    /// get data from output data
    /// First argument is file istream, and the second argument is current line istringstream
    virtual void getFromOutput(std::istream& is, std::istringstream& iss) = 0;
};

class OutputStructFilament: public OutputStruct {
public:
    static constexpr char name[] = "FILAMENT";

    OutputStructFilament(Filament* filament): _filament(filament) {}

    virtual void getFromSystemWithoutChildren()override;
    virtual void outputFromSystem(std::ostream& os)override;
    virtual void outputFromStoredWithoutChildren(std::ostream& os)override;
    virtual void getFromOutput(std::istream& is, std::istringstream& iss)override;

    int getId()const { return _id; }
    int getNumBeads()const { return _numBeads; }
    const auto& getCoords()const { return _coords; }

private:
    /// Data
    int _id;
    int _type;
    int _numBeads;
    short _deltaMinusEnd;
    short _deltaPlusEnd;

    std::vector<mathfunc::Vec< 3, floatingpoint >> _coords;

    /// Non data
    Filament* _filament;

};

class OutputStructLinker: public OutputStruct {
public:
    static constexpr char name[] = "LINKER";

    OutputStructLinker(Linker* linker): _linker(linker) {}

    virtual void getFromSystemWithoutChildren()override;
    virtual void outputFromSystem(std::ostream& os)override;
    virtual void outputFromStoredWithoutChildren(std::ostream& os)override;
    virtual void getFromOutput(std::istream& is, std::istringstream& iss)override;

    const auto& getCoords()const { return _coords; }

private:

    /// Data
    int _id;
    int _type;

    std::array<mathfunc::Vec< 3, floatingpoint >, 2> _coords;

    /// Non data
    Linker* _linker;

};

class OutputStructMotor: public OutputStruct {
public:
    static constexpr char name[] = "MOTOR";

    OutputStructMotor(MotorGhost* motor): _motor(motor) {}

    virtual void getFromSystemWithoutChildren()override;
    virtual void outputFromSystem(std::ostream& os)override;
    virtual void outputFromStoredWithoutChildren(std::ostream& os)override;
    virtual void getFromOutput(std::istream& is, std::istringstream& iss)override;

    const auto& getCoords()const { return _coords; }

private:

    /// Data
    int _id;
    int _type;
    int _bound = 1;

    std::array<mathfunc::Vec< 3, floatingpoint >, 2> _coords;

    /// Non data
    MotorGhost* _motor;

};

class OutputStructBrancher: public OutputStruct {
public:
    static constexpr char name[] = "BRANCHER";

    OutputStructBrancher(BranchingPoint* brancher): _brancher(brancher) {}

    virtual void getFromSystemWithoutChildren()override;
    virtual void outputFromSystem(std::ostream& os)override;
    virtual void outputFromStoredWithoutChildren(std::ostream& os)override;
    virtual void getFromOutput(std::istream& is, std::istringstream& iss)override;

private:

    /// Data
    int _id;
    int _type;

    mathfunc::Vec< 3, floatingpoint > _coord;

    /// Non data
    BranchingPoint* _brancher;

};

class OutputStructBubble: public OutputStruct {
public:
    static constexpr char name[] = "BUBBLE";

    OutputStructBubble(Bubble* bubble): _bubble(bubble) {}

    virtual void getFromSystemWithoutChildren()override;
    virtual void outputFromSystem(std::ostream& os)override;
    virtual void outputFromStoredWithoutChildren(std::ostream& os)override;
    virtual void getFromOutput(std::istream& is, std::istringstream& iss)override;

private:

    /// Data
    int _id;
    int _type;

    mathfunc::Vec< 3, floatingpoint > _coord;

    /// Non data
    Bubble* _bubble;

};

class OutputStructMembrane: public OutputStruct {
public:
    using Initializer = Membrane::MeshType::VertexTriangleInitializer;
    using MembraneInfo = Initializer::Info;

    static constexpr char name[] = "MEMBRANE";

    OutputStructMembrane(Membrane* membrane): _membrane(membrane) {}

    virtual void getFromSystemWithoutChildren()override;
    virtual void outputFromSystem(std::ostream& os)override;
    virtual void outputFromStoredWithoutChildren(std::ostream& os)override;
    virtual void getFromOutput(std::istream& is, std::istringstream& iss)override;

    size_t getNumVertices()const { return _numVertices; }
    size_t getNumTriangles()const { return _numTriangles; }
    const MembraneInfo& getMembraneInfo()const { return memInfo_; }

    Membrane* getMembrane()const { return _membrane; }

private:

    /// Data
    int _id;
    int _type;
    size_t _numVertices;
    size_t _numTriangles;
    size_t numBorders_ = 0;

    MembraneInfo memInfo_;

    /// Non data
    Membrane* _membrane;
};

class OutputStructSnapshot: public OutputStruct {
private:
    /// Bonus data (generated when constructed)
    int _snapshot;

    /// Data
    int _numFilaments;
    int _numLinkers;
    int _numMotorGhosts;
    int _numBranchingPoints;
    int _numBubbles;
    int _numMembranes;

public:
    // data
    double simulationTime = 0;

    /// Children data
    std::vector<OutputStructFilament> filamentStruct;
    std::vector<OutputStructLinker>   linkerStruct;
    std::vector<OutputStructMotor>    motorStruct;
    std::vector<OutputStructBrancher> brancherStruct;
    std::vector<OutputStructBubble>   bubbleStruct;
    std::vector<OutputStructMembrane> membraneStruct;

    /// Non data

    OutputStructSnapshot(int snapshot): _snapshot(snapshot) {}

    virtual void getFromSystem()override;
    virtual void getFromSystemWithoutChildren()override;
    virtual void outputFromSystem(std::ostream& os)override;
    virtual void outputFromStored(std::ostream& os)override;
    virtual void outputFromStoredWithoutChildren(std::ostream& os)override;
    virtual void getFromOutput(std::istream& is, std::istringstream& iss)override;
};

#endif
