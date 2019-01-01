#include "Structure/SurfaceMesh/Membrane.hpp"

#include <unordered_set>

#include "Compartment.h"
#include "core/controller/GController.h"
#include "SubSystem.h"
#include "SysParams.h"

#include "Triangle.h"
#include "Edge.h"
#include "Vertex.h"
#include "MTriangle.h"
#include "MVoronoiCell.h"
#include "MMembrane.h"

using namespace mathfunc;

Database<Membrane*> Membrane::_membranes;

Membrane::Membrane(
    SubSystem* s,
    short membraneType,
    const std::vector< MembraneMeshAttribute::coordinate_type >& vertexCoordinateList,
    const std::vector< std::array< size_t, 3 > >& triangleVertexIndexList
) : Trackable(false, false, false, false),
    _mesh(MembraneMeshAttribute::MetaAttribute{s, this}),
    _subSystem(s), _memType(membraneType), _id(_membranes.getID()) {
    
    // Build the meshwork using vertex and triangle information
    _mesh.init(vertexCoordinateList, triangleVertexIndexList);

    // Update geometry
    updateGeometryValue<false>();

    size_t numTriangles = _mesh.getTriangles().size();

#ifdef MECHANICS
    // Calculate the total area and volume to set the equilibrium area and volume
    double area = 0.0;
    double volume = 0.0;
    for(const auto& t : _mesh.getTriangles()) {
        area += t.attr.gTriangle.area;
        volume += t.attr.gTriangle.coneVolume;
    }

    _mMembrane = std::make_unique<MMembrane>(_memType);
    _mMembrane->setEqArea(
        area * SysParams::Mechanics().MemEqAreaFactor[membraneType]
    );
    _mMembrane->setEqVolume(volume);
#endif

}

void Membrane::printSelf()const {
    
    cout << endl;
    
    cout << "Membrane: ptr = " << this << endl;
    cout << "Membrane Id = " << _id << endl;
    cout << "Membrane type = " << _memType << endl;
    
    cout << endl;
    
}

double Membrane::meshworkQuality()const {
    /*
    This function calculates the quality of the meshwork of this membrane, and
    the result is represented as a value between 0 and 1, 0 being worst and 1
    being the best case (all equilateral).

    The criterion used is the minimum angle in each triangle, parametrized and
    averaged over all triangles. And the calculation requires the result of
        - The angle calculation of all triangles
    */

    double res = 0;

    // TODO implementation

    return res;
}
