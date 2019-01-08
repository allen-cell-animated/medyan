
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
#include "OutputStruct.h"

#include <iterator>
#include <string>

#include "common.h"
#include "MathFunctions.h"

#include "Parser.h"

#include "Bead.h"
#include "BranchingPoint.h"
#include "Bubble.h"
#include "Cylinder.h"
#include "Filament.h"
#include "Linker.h"
#include "Membrane.hpp"
#include "MotorGhost.h"
#include "Vertex.h"

using namespace mathfunc;

/******************************************************************************
OutputStruct for Filaments
******************************************************************************/
//@{
constexpr char OutputStructFilament::name[];

void OutputStructFilament::getFromSystemWithoutChildren() {
    _id = _filament->getID();
    _type = _filament->getType();
    _numBeads = _filament->getCylinderVector().size() + 1;
    _deltaMinusEnd = _filament->getDeltaMinusEnd();
    _deltaPlusEnd = _filament->getDeltaPlusEnd();

    // Store coordinates
    _coords.reserve(_numBeads);
    for (auto cylinder : _filament->getCylinderVector())
        _coords.push_back(vector2Vec<3, double>(cylinder->getFirstBead()->coordinate));
    _coords.push_back(vector2Vec<3, double>(_filament->getCylinderVector().back()->getSecondBead()->coordinate));
        
}

void OutputStructFilament::outputFromSystem(std::ostream& os) {
    getFromSystemWithoutChildren();
    outputFromStoredWithoutChildren(os);
}

void OutputStructFilament::outputFromStoredWithoutChildren(std::ostream& os) {
    os << name << " "
        << _id << " "
        << _type << " "
        << _numBeads << " "
        << _deltaMinusEnd << " "
        << _deltaPlusEnd << std::endl;

    for(auto& coord: _coords)
        for(double value: coord)
            os << value << " ";
    os << std::endl;
}

void OutputStructFilament::getFromOutput(std::istream& is, std::istringstream& iss) {
    iss >> _id
        >> _type
        >> _numBeads
        >> _deltaMinusEnd
        >> _deltaPlusEnd;

    std::string nextLine;
    std::getline(is, nextLine);
    std::istringstream newIss(nextLine);
    _coords.clear();
    double tmp;
    while(newIss >> tmp) {
        mathfunc::Vec3 coord;
        coord[0] = tmp;
        newIss >> coord[1] >> coord[2];
        _coords.push_back(coord);
    }
}
//@}

/******************************************************************************
OutputStruct for Linkers
******************************************************************************/
//@{
constexpr char OutputStructLinker::name[];

void OutputStructLinker::getFromSystemWithoutChildren() {
    _id = _linker->getID();
    _type = _linker->getType();

    // Store coordinates
    _coords = {{
        vector2Vec<3, double>(midPointCoordinate(
            _linker->getFirstCylinder()->getFirstBead()->coordinate,
            _linker->getFirstCylinder()->getSecondBead()->coordinate,
            _linker->getFirstPosition()
        )),
        vector2Vec<3, double>(midPointCoordinate(
            _linker->getSecondCylinder()->getFirstBead()->coordinate,
            _linker->getSecondCylinder()->getSecondBead()->coordinate,
            _linker->getSecondPosition()
        ))
    }};
        
}

void OutputStructLinker::outputFromSystem(std::ostream& os) {
    getFromSystemWithoutChildren();
    outputFromStoredWithoutChildren(os);
}

void OutputStructLinker::outputFromStoredWithoutChildren(std::ostream& os) {
    os << name << " "
        << _id << " "
        << _type << std::endl;

    for(auto& coord: _coords)
        for(double value: coord)
            os << value << " ";
    os << std::endl;
}

void OutputStructLinker::getFromOutput(std::istream& is, std::istringstream& iss) {
    iss >> _id
        >> _type;

    std::string nextLine;
    std::getline(is, nextLine);
    std::istringstream newIss(nextLine);
    for(auto& coord: _coords)
        for(double& value: coord)
            newIss >> value;
}
//@}

/******************************************************************************
OutputStruct for Motors
******************************************************************************/
//@{
constexpr char OutputStructMotor::name[];

void OutputStructMotor::getFromSystemWithoutChildren() {
    _id = _motor->getID();
    _type = _motor->getType();

    // Store coordinates
    _coords = {{
        vector2Vec<3, double>(midPointCoordinate(
            _motor->getFirstCylinder()->getFirstBead()->coordinate,
            _motor->getFirstCylinder()->getSecondBead()->coordinate,
            _motor->getFirstPosition()
        )),
        vector2Vec<3, double>(midPointCoordinate(
            _motor->getSecondCylinder()->getFirstBead()->coordinate,
            _motor->getSecondCylinder()->getSecondBead()->coordinate,
            _motor->getSecondPosition()
        ))
    }};
        
}

void OutputStructMotor::outputFromSystem(std::ostream& os) {
    getFromSystemWithoutChildren();
    outputFromStoredWithoutChildren(os);
}

void OutputStructMotor::outputFromStoredWithoutChildren(std::ostream& os) {
    os << name << " "
        << _id << " "
        << _type << " "
        << _bound << std::endl;

    for(auto& coord: _coords)
        for(double value: coord)
            os << value << " ";
    os << std::endl;
}

void OutputStructMotor::getFromOutput(std::istream& is, std::istringstream& iss) {
    iss >> _id
        >> _type
        >> _bound;

    std::string nextLine;
    std::getline(is, nextLine);
    std::istringstream newIss(nextLine);
    for(auto& coord: _coords)
        for(double& value: coord)
            newIss >> value;
}
//@}

/******************************************************************************
OutputStruct for Branchers
******************************************************************************/
//@{
constexpr char OutputStructBrancher::name[];

void OutputStructBrancher::getFromSystemWithoutChildren() {
    _id = _brancher->getID();
    _type = _brancher->getType();

    // Store coordinates
    _coord = vector2Vec<3, double>(_brancher->coordinate);
        
}

void OutputStructBrancher::outputFromSystem(std::ostream& os) {
    getFromSystemWithoutChildren();
    outputFromStoredWithoutChildren(os);
}

void OutputStructBrancher::outputFromStoredWithoutChildren(std::ostream& os) {
    os << name << " "
        << _id << " "
        << _type << std::endl;

    for(double value: _coord)
        os << value << " ";
    os << std::endl;
}

void OutputStructBrancher::getFromOutput(std::istream& is, std::istringstream& iss) {
    iss >> _id
        >> _type;

    std::string nextLine;
    std::getline(is, nextLine);
    std::istringstream newIss(nextLine);
    for(double& value: _coord)
        newIss >> value;
}
//@}

/******************************************************************************
OutputStruct for Bubbles
******************************************************************************/
//@{
constexpr char OutputStructBubble::name[];

void OutputStructBubble::getFromSystemWithoutChildren() {
    _id = _bubble->getID();
    _type = _bubble->getType();

    // Store coordinates
    _coord = vector2Vec<3, double>(_bubble->coordinate);
        
}

void OutputStructBubble::outputFromSystem(std::ostream& os) {
    getFromSystemWithoutChildren();
    outputFromStoredWithoutChildren(os);
}

void OutputStructBubble::outputFromStoredWithoutChildren(std::ostream& os) {
    os << name << " "
        << _id << " "
        << _type << std::endl;

    for(double value: _coord)
        os << value << " ";
    os << std::endl;
}

void OutputStructBubble::getFromOutput(std::istream& is, std::istringstream& iss) {
    iss >> _id
        >> _type;

    std::string nextLine;
    std::getline(is, nextLine);
    std::istringstream newIss(nextLine);
    for(double& value: _coord)
        newIss >> value;
}
//@}

/******************************************************************************
OutputStruct for Membranes
******************************************************************************/
//@{
constexpr char OutputStructMembrane::name[];

void OutputStructMembrane::getFromSystemWithoutChildren() {
    _id = _membrane->getId();
    _type = _membrane->getType();

    _memInfo = _membrane->getMesh().extract< Initializer >();
    _numVertices = _memInfo.attributeInitializerInfo.vertexCoordinateList.size();
    _numTriangles = _memInfo.triangleVertexIndexList.size();
        
}

void OutputStructMembrane::outputFromSystem(std::ostream& os) {
    getFromSystemWithoutChildren();
    outputFromStoredWithoutChildren(os);
}

void OutputStructMembrane::outputFromStoredWithoutChildren(std::ostream& os) {
    os << name << " "
        << _id << " "
        << _type << " "
        << _numVertices << ' '
        << _numTriangles << '\n';

    // print coordinates
    for(const auto& it : _memInfo.attributeInitializerInfo.vertexCoordinateList) {
        for(double value : it) os << value << ' ';
        os << '\n';
    }

    // print neighbor indices
    for(const auto& it : _memInfo.triangleVertexIndexList) {
        for(size_t value : it) os << value << ' ';
        os << '\n';
    }
}

void OutputStructMembrane::getFromOutput(std::istream& is, std::istringstream& iss) {
    iss >> _id
        >> _type
        >> _numVertices
        >> _numTriangles;

    _memInfo.attributeInitializerInfo.vertexCoordinateList.reserve(_numVertices);
    _memInfo.triangleVertexIndexList.reserve(_numTriangles);
    for(size_t i = 0; i < _numVertices; ++i) {
        std::string nextLine;
        std::getline(is, nextLine);
        std::istringstream newIss(nextLine);

        _memInfo.attributeInitializerInfo.vertexCoordinateList.emplace_back(3);
        auto& coord = _memInfo.attributeInitializerInfo.vertexCoordinateList.back();

        // Record coordinates
        for(double& value: coord)
            newIss >> value;
    }
    for(size_t i = 0; i < _numTriangles; ++i) {
        std::string nextLine;
        std::getline(is, nextLine);
        std::istringstream newIss(nextLine);

        _memInfo.triangleVertexIndexList.emplace_back();
        auto& indices = _memInfo.triangleVertexIndexList.back();

        // Record indices
        for(size_t& value: indices)
            newIss >> value;
    }
}

//@}

/******************************************************************************
OutputStruct for snapshots
******************************************************************************/
//@{
void OutputStructSnapshot::getFromSystem() {
    getFromSystemWithoutChildren();

    for(auto filament: Filament::getFilaments()) {
        filamentStruct.emplace_back(filament);
        filamentStruct.back().getFromSystem();
    }
    for(auto linker: Linker::getLinkers()) {
        linkerStruct.emplace_back(linker);
        linkerStruct.back().getFromSystem();
    }
    for(auto motor: MotorGhost::getMotorGhosts()) {
        motorStruct.emplace_back(motor);
        motorStruct.back().getFromSystem();
    }
    for(auto brancher: BranchingPoint::getBranchingPoints()) {
        brancherStruct.emplace_back(brancher);
        brancherStruct.back().getFromSystem();
    }
    for(auto bubble: Bubble::getBubbles()) {
        bubbleStruct.emplace_back(bubble);
        bubbleStruct.back().getFromSystem();
    }
    for(auto membrane: Membrane::getMembranes()) {
        membraneStruct.emplace_back(membrane);
        membraneStruct.back().getFromSystem();
    }

    // Note: new children should be added here
}

void OutputStructSnapshot::getFromSystemWithoutChildren() {
    _time = tau();

    _numFilaments = Filament::numFilaments();
    _numLinkers = Linker::numLinkers();
    _numMotorGhosts = MotorGhost::numMotorGhosts();
    _numBranchingPoints = BranchingPoint::numBranchingPoints();
    _numBubbles = Bubble::numBubbles();
    _numMembranes = Membrane::numMembranes();

}

void OutputStructSnapshot::outputFromSystem(std::ostream& os) {
    getFromSystemWithoutChildren();
    outputFromStoredWithoutChildren(os);

    for(auto filament: Filament::getFilaments()) {
        filamentStruct.emplace_back(filament);
        filamentStruct.back().outputFromSystem(os);
    }
    for(auto linker: Linker::getLinkers()) {
        linkerStruct.emplace_back(linker);
        linkerStruct.back().outputFromSystem(os);
    }
    for(auto motor: MotorGhost::getMotorGhosts()) {
        motorStruct.emplace_back(motor);
        motorStruct.back().outputFromSystem(os);
    }
    for(auto brancher: BranchingPoint::getBranchingPoints()) {
        brancherStruct.emplace_back(brancher);
        brancherStruct.back().outputFromSystem(os);
    }
    for(auto bubble: Bubble::getBubbles()) {
        bubbleStruct.emplace_back(bubble);
        bubbleStruct.back().outputFromSystem(os);
    }
    for(auto membrane: Membrane::getMembranes()) {
        membraneStruct.emplace_back(membrane);
        membraneStruct.back().outputFromSystem(os);
    }

    // Note: new children should be added here
}

void OutputStructSnapshot::outputFromStored(std::ostream& os) {
    outputFromStoredWithoutChildren(os);
    
    for(auto& it: filamentStruct) {
        it.outputFromStored(os);
    }
    for(auto& it: linkerStruct) {
        it.outputFromStored(os);
    }
    for(auto& it: motorStruct) {
        it.outputFromStored(os);
    }
    for(auto& it: brancherStruct) {
        it.outputFromStored(os);
    }
    for(auto& it: bubbleStruct) {
        it.outputFromStored(os);
    }
    for(auto& it: membraneStruct) {
        it.outputFromStored(os);
    }

    // Note: new children should be added here
}

void OutputStructSnapshot::outputFromStoredWithoutChildren(std::ostream& os) {
    os << _snapshot << " "
        << _time << " "
        << _numFilaments << " "
        << _numLinkers << " "
        << _numMotorGhosts << " "
        << _numBranchingPoints << " "
        << _numBubbles << " "
        << _numMembranes << std::endl;
}

void OutputStructSnapshot::getFromOutput(std::istream& is, std::istringstream& iss) {
    iss >> _snapshot
        >> _time
        >> _numFilaments
        >> _numLinkers
        >> _numMotorGhosts
        >> _numBranchingPoints
        >> _numBubbles
        >> _numMembranes;

    std::string nextLine;
    do {
        std::getline(is, nextLine);
        if(!is || nextLine.empty()) break;

        std::istringstream newIss(nextLine);
        std::string name;
        newIss >> name;
        if(name == OutputStructFilament::name) {
            filamentStruct.emplace_back(nullptr);
            filamentStruct.back().getFromOutput(is, newIss);
        } else if(name == OutputStructLinker::name) {
            linkerStruct.emplace_back(nullptr);
            linkerStruct.back().getFromOutput(is, newIss);
        } else if(name == OutputStructMotor::name) {
            motorStruct.emplace_back(nullptr);
            motorStruct.back().getFromOutput(is, newIss);
        } else if(name == OutputStructBrancher::name) {
            brancherStruct.emplace_back(nullptr);
            brancherStruct.back().getFromOutput(is, newIss);
        } else if(name == OutputStructBubble::name) {
            bubbleStruct.emplace_back(nullptr);
            bubbleStruct.back().getFromOutput(is, newIss);
        } else if(name == OutputStructMembrane::name) {
            membraneStruct.emplace_back(nullptr);
            membraneStruct.back().getFromOutput(is, newIss);
        } // TODO: other children

    } while(true);
}
//@}
