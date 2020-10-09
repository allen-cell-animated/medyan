#include "VisualSystemRawData.hpp"

#include "Structure/BranchingPoint.h"
#include "Structure/Filament.h"
#include "Structure/Linker.h"
#include "Structure/MotorGhost.h"
#include "Structure/SurfaceMesh/Membrane.hpp"

namespace medyan::visual {

bool copySystemData(
    SystemRawData& data,
    raw_data_cat::Type updated,
    bool ignoreDataInUse
) {
    std::unique_lock lk(data.me, std::try_to_lock);
    if(!lk.owns_lock()) {
        if(ignoreDataInUse) {
            // The data is in use, skip copying
            return false;
        } else {
            // Acquire the lock
            lk.lock();
        }
    }

    // update simulation time
    data.frameData.simulationTime = tau();

    if(updated & (raw_data_cat::beadPosition | raw_data_cat::beadConnection)) {

        // Extract membrane indexing
        data.frameData.membranes.clear();
        for(const Membrane* m : Membrane::getMembranes()) {
            const auto& mesh = m->getMesh();

            data.frameData.membranes.emplace_back();
            auto& mi = data.frameData.membranes.back();

            mi.vertexCoords.reserve(mesh.numVertices());
            for(const auto& v : mesh.getVertices()) {
                mi.vertexCoords.push_back(mathfunc::Vec3f(v.attr.vertex->coord));
            }

            mi.triangles.reserve(mesh.numTriangles());
            for(const auto& t : mesh.getTriangles()) {
                size_t vIndex = 0; // 0, 1, 2
                mi.triangles.emplace_back();
                mesh.forEachHalfEdgeInPolygon(t, [&](auto hei) {
                    mi.triangles.back()[vIndex] = mesh.target(hei).index;
                    ++vIndex;
                });
            }
        }

        // Extract filament indexing
        data.frameData.filaments.clear();
        for(Filament* f : Filament::getFilaments()) {
            data.frameData.filaments.emplace_back();
            auto& fi = data.frameData.filaments.back().coords;

            const auto& cylinders = f->getCylinderVector();
            fi.reserve(cylinders.size() + 1);
            for(Cylinder* c : cylinders)
                fi.push_back(Vec3f(c->getFirstBead()->coordinate()));
            fi.push_back(Vec3f(cylinders.back()->getSecondBead()->coordinate()));
        }

        // Extract motors, linkers and branchers
        data.frameData.linkers.clear();
        for(Linker* l : Linker::getLinkers()) {
            auto& newLinker = data.frameData.linkers.emplace_back();

            const auto pos0 = l->getFirstPosition();
            newLinker.coords[0]
                = l->getFirstCylinder()->getFirstBead()->coordinate() * (1 - pos0)
                + l->getFirstCylinder()->getSecondBead()->coordinate() * pos0;
            const auto pos1 = l->getSecondPosition();
            newLinker.coords[1]
                = l->getSecondCylinder()->getFirstBead()->coordinate() * (1 - pos1)
                + l->getSecondCylinder()->getSecondBead()->coordinate() * pos1;
            
            newLinker.type = displayLinkerType(data.displayTypeMap, "linker");
        }

        for(MotorGhost* l : MotorGhost::getMotorGhosts()) {
            auto& newLinker = data.frameData.linkers.emplace_back();

            const auto pos0 = l->getFirstPosition();
            newLinker.coords[0]
                = l->getFirstCylinder()->getFirstBead()->coordinate() * (1 - pos0)
                + l->getFirstCylinder()->getSecondBead()->coordinate() * pos0;
            const auto pos1 = l->getSecondPosition();
            newLinker.coords[1]
                = l->getSecondCylinder()->getFirstBead()->coordinate() * (1 - pos1)
                + l->getSecondCylinder()->getSecondBead()->coordinate() * pos1;

            newLinker.type = displayLinkerType(data.displayTypeMap, "motor");
        }

        for(BranchingPoint* l : BranchingPoint::getBranchingPoints()) {
            auto& newLinker = data.frameData.linkers.emplace_back();

            const auto pos = l->getPosition();
            newLinker.coords[0]
                = l->getFirstCylinder()->getFirstBead()->coordinate() * (1 - pos)
                + l->getFirstCylinder()->getSecondBead()->coordinate() * pos;
            newLinker.coords[1]
                = l->getSecondCylinder()->getFirstBead()->coordinate();

            newLinker.type = displayLinkerType(data.displayTypeMap, "brancher");
        }

    } // End update bead connection

    if(updated & (raw_data_cat::compartment)) {
        auto& compartmentInfo = data.frameData.compartmentInfo;
        compartmentInfo.emplace();

        compartmentInfo->number = {
            SysParams::Geometry().NX,
            SysParams::Geometry().NY,
            SysParams::Geometry().NZ
        };
        compartmentInfo->size = {
            static_cast< float >(SysParams::Geometry().compartmentSizeX),
            static_cast< float >(SysParams::Geometry().compartmentSizeY),
            static_cast< float >(SysParams::Geometry().compartmentSizeZ)
        };

    }

    // Save updated
    data.updated |= updated;

    return true;
}

} // namespace medyan::visual
