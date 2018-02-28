
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
#include "Membrane.h"
#include "MotorGhost.h"
#include "Vertex.h"

using namespace mathfunc;

/******************************************************************************
OutputStruct for Filaments
******************************************************************************/
//@{
void OutputStructFilament::getFromSystemWithoutChildren() {
    _id = _filament->getID();
    _type = _filament->getType();
    _numBeads = _filament->getCylinderVector().size() + 1;
    _deltaMinusEnd = _filament->getDeltaMinusEnd();
    _deltaPlusEnd = _filament->getDeltaPlusEnd();

    // Store coordinates
    _coords.reserve(_numBeads);
    for (auto cylinder : _filament->getCylinderVector())
        _coords.push_back(vector2Array<double, 3>(cylinder->getFirstBead()->coordinate));
    _coords.push_back(vector2Array<double, 3>(_filament->getCylinderVector().back()->getSecondBead()->coordinate));
        
    // Side effect: Reset deltas for this filament
    // TODO: Move to output
    _filament->resetDeltaPlusEnd();
    _filament->resetDeltaMinusEnd();
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
        std::array<double, 3> coord;
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
void OutputStructLinker::getFromSystemWithoutChildren() {
    _id = _linker->getID();
    _type = _linker->getType();

    // Store coordinates
    _coords = {{
        vector2Array<double, 3>(midPointCoordinate(
            _linker->getFirstCylinder()->getFirstBead()->coordinate,
            _linker->getFirstCylinder()->getSecondBead()->coordinate,
            _linker->getFirstPosition()
        )),
        vector2Array<double, 3>(midPointCoordinate(
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
void OutputStructMotor::getFromSystemWithoutChildren() {
    _id = _motor->getID();
    _type = _motor->getType();

    // Store coordinates
    _coords = {{
        vector2Array<double, 3>(midPointCoordinate(
            _motor->getFirstCylinder()->getFirstBead()->coordinate,
            _motor->getFirstCylinder()->getSecondBead()->coordinate,
            _motor->getFirstPosition()
        )),
        vector2Array<double, 3>(midPointCoordinate(
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
void OutputStructBrancher::getFromSystemWithoutChildren() {
    _id = _brancher->getID();
    _type = _brancher->getType();

    // Store coordinates
    _coord = vector2Array<double, 3>(_brancher->coordinate);
        
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
void OutputStructBubble::getFromSystemWithoutChildren() {
    _id = _bubble->getID();
    _type = _bubble->getType();

    // Store coordinates
    _coord = vector2Array<double, 3>(_bubble->coordinate);
        
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
void OutputStructMembrane::getFromSystemWithoutChildren() {
    _id = _membrane->getId();
    _type = _membrane->getType();

    auto& vertices = _membrane->getVertexVector();
    _numVertices = vertices.size();

    // Store coordinates with neighbor indices
    _memInfo.reserve(_numVertices);
    for(size_t idx = 0; idx < _numVertices; ++idx) {

        _memInfo.emplace_back();
        auto& vtxInfo = _memInfo.back();
        auto& coord = get<0>(vtxInfo);
        auto& neighborIndices = get<1>(vtxInfo);
            
        coord = vector2Array<double, 3>(vertices[idx]->coordinate);
        
        auto& neighbors = vertices[idx]->getNeighborVertices();
        size_t numNeighbors = neighbors.size();
        neighborIndices.reserve(numNeighbors);

        for(size_t nIdx = 0; nIdx < numNeighbors; ++nIdx) {
            neighborIndices.push_back(neighbors[nIdx]->getMembraneVertexIdx());
        }
    }
        
}

void OutputStructMembrane::outputFromSystem(std::ostream& os) {
    getFromSystemWithoutChildren();
    outputFromStoredWithoutChildren(os);
}

void OutputStructMembrane::outputFromStoredWithoutChildren(std::ostream& os) {
    os << name << " "
        << _id << " "
        << _type << " "
        << _numVertices << std::endl;

    for(auto& vtxInfo: _memInfo) {
        
        // print coordinates
        for(double value: get<0>(vtxInfo))
            os << value << " ";

        // print neighbor indices
        for(size_t value: get<1>(vtxInfo))
            os << value << " ";

        os << std::endl;
        
    }
}

void OutputStructMembrane::getFromOutput(std::istream& is, std::istringstream& iss) {
    iss >> _id
        >> _type
        >> _numVertices;

    _memInfo.reserve(_numVertices);
    for(size_t idx = 0; idx < _numVertices; ++idx) {
        std::string nextLine;
        std::getline(is, nextLine);
        std::istringstream newIss(nextLine);

        _memInfo.emplace_back();
        auto& vtxInfo = _memInfo.back();
        auto& coord = get<0>(vtxInfo);
        auto& neighborIndices = get<1>(vtxInfo);

        // Record coordinates
        for(double& value: coord)
            newIss >> value;
        
        // Record neighbor indices
        neighborIndices = vector<size_t>(istream_iterator<size_t>(newIss), istream_iterator<size_t>());

    }
}

size_t OutputStructMembrane::getNumEdges()const {
    if(_membrane) {
        return _membrane->getEdgeVector().size();
    } else {
        size_t numEdges = 0;
        for(const VertexInfo& v: _memInfo) {
            numEdges += get<1>(v).size(); // Add number of neighbors
        }
        numEdges /= 2;
        return numEdges;
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
