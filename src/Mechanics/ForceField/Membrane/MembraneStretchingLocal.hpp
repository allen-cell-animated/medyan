#ifndef MEDYAN_Mechanics_ForceField_Membrane_MembraneStretching_hpp
#define MEDYAN_Mechanics_ForceField_Membrane_MembraneStretching_hpp

#include <array>
#include <limits>
#include <type_traits>
#include <vector>

#include "Mechanics/ForceField/ForceField.h"
#include "Mechanics/ForceField/Membrane/MembraneStretchingImpl.hpp"
#include "Structure/SurfaceMesh/Membrane.hpp"

enum class MembraneStretchingType {
    localHarmonic,
    globalHarmonic
};

struct MembraneStretchingLocal : public ForceField {

    MembraneStretchingHarmonic impl;

    // (temp) holds the first index of coordinates of each vertex of triangles
    std::vector< std::array< std::size_t, 3 >> vertexSet;

    // (temp) holds the mech params of triangles
    std::vector< double > kArea;
    std::vector< double > eqArea;

    virtual void vectorize(const FFCoordinateStartingIndex& si) override {
        using namespace std;
        using AT = Membrane::MembraneMeshAttributeType;

        vertexSet.clear();
        vertexSet.reserve(Triangle::numElements()); // Might be more than needed
        kArea.clear();
        kArea.reserve(Triangle::numElements());
        eqArea.clear();
        eqArea.reserve(Triangle::numElements());

        for(auto m : Membrane::getMembranes()) {
            // In area elasticity for each triangle, the triangles with any
            // vertex touching an open border will not be included.
            // TODO: add triangle filter

            const auto& mesh = m->getMesh();
            for(const auto& t : mesh.getTriangles()) {
                // TODO: check if the triangle fits local/global criteria
                const auto vis = AT::vertexIndices(mesh, t);
                vertexSet.push_back({
                    mesh.attribute(vis[0]).vertex->getIndex() * 3 + si.vertex,
                    mesh.attribute(vis[1]).vertex->getIndex() * 3 + si.vertex,
                    mesh.attribute(vis[2]).vertex->getIndex() * 3 + si.vertex
                });
                kArea.push_back(t.attr.triangle->mTriangle.kArea);
                eqArea.push_back(t.attr.triangle->mTriangle.eqArea);
            }
        }
    }

    virtual floatingpoint computeEnergy(floatingpoint* coord, bool stretched) override {
        using namespace std;
        using AT = Membrane::MembraneMeshAttributeType;

        double en = 0;

        for(size_t i = 0; i < vertexSet.size(); ++i) {
            const auto& vs = vertexSet[i];
            const auto enTriangle = impl.energy(
                AT::area(
                    makeRefVec<3>(coord + vs[0]),
                    makeRefVec<3>(coord + vs[1]),
                    makeRefVec<3>(coord + vs[2])
                ),
                kArea[i],
                eqArea[i]
            );

            if(!isfinite(enTriangle)) {
                LOG(ERROR) << "In " << getName() << " energy calculation, triangle with index "
                    << i << " has energy " << enTriangle;

                return numeric_limits<double>::infinity();
            }
            else {
                en += enMem;
            }

        }

        return en;
    }

    virtual void computeForces(floatingpoint* coord, floatingpoint* force) override {
        using namespace std;
        using AT = Membrane::MembraneMeshAttributeType;

        for(size_t i = 0; i < vertexSet.size(); ++i) {
            const auto& vs = vertexSet[i];
            const auto rv0 = makeRefVec<3>(coord + vs[0]);
            const auto rv1 = makeRefVec<3>(coord + vs[1]);
            const auto rv2 = makeRefVec<3>(coord + vs[2]);
            const auto [area, da0, da1, da2] = AT::areaAndDerivative(rv0, rv1, rv2);

            impl.forces(force + vs[0], area, da0, kArea[i], eqArea[i]);
            impl.forces(force + vs[1], area, da1, kArea[i], eqArea[i]);
            impl.forces(force + vs[2], area, da2, kArea[i], eqArea[i]);

        }

    }

    virtual string getName() override { return "Membrane Stretching local"; }


    // Useless overrides
    virtual void cleanup() override {}
    virtual void computeLoadForces() override {}
    virtual void whoIsCulprit() override {}
    virtual std::vector<NeighborList*> getNeighborLists() override { return {}; }
};

#endif
