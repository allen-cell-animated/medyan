#include "Structure/SurfaceMesh/Membrane.hpp"

#include "Compartment.h"
#include "core/controller/GController.h"
#include "SubSystem.h"
#include "SysParams.h"

#include "MTriangle.h"
#include "MVoronoiCell.h"
#include "MMembrane.h"

using namespace mathfunc;

Database<Membrane*> Membrane::_membranes;

Membrane::Membrane(
    SubSystem* s,
    short membraneType,
    const std::vector< coordinate_type >& vertexCoordinateList,
    const std::vector< std::array< size_t, 3 > >& triangleVertexIndexList
) : Trackable(false, false, false, false),
    _mesh(typename MembraneMeshAttributeType::MetaAttribute{s, this}),
    _subSystem(s), _memType(membraneType), _id(_membranes.getID()) {
    
    // Build the meshwork topology using vertex and triangle information
    _mesh.init<typename MeshType::VertexTriangleInitializer>(
        vertexCoordinateList.size(),
        triangleVertexIndexList,
        typename MembraneMeshAttributeType::AttributeInitializerInfo{ vertexCoordinateList }
    );

    // Update geometry
    updateGeometryValue<false>();

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
