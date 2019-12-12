#include "VisualSystemRawData.hpp"

#include "Structure/BranchingPoint.h"
#include "Structure/Filament.h"
#include "Structure/Linker.h"
#include "Structure/MotorGhost.h"
#include "Structure/SurfaceMesh/Membrane.hpp"

namespace visual {

void copySystemData(
    SystemRawData& data,
    raw_data_cat::Type updated,
    bool ignoreDataInUse
) {
    std::unique_lock lk(data.me, std::try_to_lock);
    if(!lk.owns_lock() && ignoreDataInUse) {
        // The data is in use, skip copying
        return;
    }

    if(updated & (raw_data_cat::beadPosition | raw_data_cat::beadConnection)) {
        // Copy bead data
        data.copiedBeadData = Bead::getDbDataConst();
    }

    if(updated & (raw_data_cat::beadConnection)) {

        // Extract membrane indexing
        data.membraneIndices.clear();
        for(const Membrane* m : Membrane::getMembranes()) {
            const auto& mesh = m->getMesh();

            data.membraneIndices.emplace_back();
            auto& mi = data.membraneIndices.back();

            mi.vertexIndices.reserve(mesh.numVertices());
            for(const auto& v : mesh.getVertices()) {
                mi.vertexIndices.push_back(v.attr.vertex->Bead::getStableIndex());
            }

            mi.triangleVertexIndices.reserve(mesh.numTriangles());
            for(const auto& t : mesh.getTriangles()) {
                size_t vIndex = 0; // 0, 1, 2
                mi.triangleVertexIndices.emplace_back();
                mesh.forEachHalfEdgeInPolygon(t, [&](size_t hei) {
                    mi.triangleVertexIndices.back()[vIndex] = mesh.target(hei);
                    ++vIndex;
                });
            }
        }

        // Extract filament indexing
        data.filamentIndices.clear();
        for(Filament* f : Filament::getFilaments()) {
            data.filamentIndices.emplace_back();
            auto& fi = data.filamentIndices.back();

            const auto& cylinders = f->getCylinderVector(); // TODO make f const Filament*
            fi.reserve(cylinders.size() + 1);
            for(Cylinder* c : cylinders)
                fi.push_back(c->getFirstBead()->getStableIndex()); // TODO make c const Cylinder*
            fi.push_back(cylinders.back()->getSecondBead()->getStableIndex());
        }

        // Extract motors, linkers and branchers
        data.linkerCoords.clear();
        for(Linker* l : Linker::getLinkers()) {
            data.linkerCoords.emplace_back();
            const auto pos0 = l->getFirstPosition();
            data.linkerCoords.back()[0]
                = l->getFirstCylinder()->getFirstBead()->coordinate() * (1 - pos0)
                + l->getFirstCylinder()->getSecondBead()->coordinate() * pos0;
            const auto pos1 = l->getSecondPosition();
            data.linkerCoords.back()[1]
                = l->getSecondCylinder()->getFirstBead()->coordinate() * (1 - pos1)
                + l->getSecondCylinder()->getSecondBead()->coordinate() * pos1;
        }
        data.motorCoords.clear();
        for(MotorGhost* l : MotorGhost::getMotorGhosts()) {
            data.motorCoords.emplace_back();
            const auto pos0 = l->getFirstPosition();
            data.motorCoords.back()[0]
                = l->getFirstCylinder()->getFirstBead()->coordinate() * (1 - pos0)
                + l->getFirstCylinder()->getSecondBead()->coordinate() * pos0;
            const auto pos1 = l->getSecondPosition();
            data.motorCoords.back()[1]
                = l->getSecondCylinder()->getFirstBead()->coordinate() * (1 - pos1)
                + l->getSecondCylinder()->getSecondBead()->coordinate() * pos1;
        }
        data.brancherCoords.clear();
        for(BranchingPoint* l : BranchingPoint::getBranchingPoints()) {
            data.brancherCoords.emplace_back();
            const auto pos = l->getPosition();
            data.brancherCoords.back()[0]
                = l->getFirstCylinder()->getFirstBead()->coordinate() * (1 - pos)
                + l->getFirstCylinder()->getSecondBead()->coordinate() * pos;
            data.brancherCoords.back()[1]
                = l->getSecondCylinder()->getFirstBead()->coordinate();
        }

    } // End update bead connection

    if(updated & (raw_data_cat::compartment)) {
        data.compartmentNum = {
            static_cast< size_t >(SysParams::Geometry().NX),
            static_cast< size_t >(SysParams::Geometry().NY),
            static_cast< size_t >(SysParams::Geometry().NZ)
        };
        data.compartmentSize = {
            SysParams::Geometry().compartmentSizeX,
            SysParams::Geometry().compartmentSizeY,
            SysParams::Geometry().compartmentSizeZ
        };

    }

    // Save updated
    data.updated |= updated;

}

} // namespace visual
