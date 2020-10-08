#ifndef MEDYAN_Mechanics_ForceField_Membrane_MembraneTension_hpp
#define MEDYAN_Mechanics_ForceField_Membrane_MembraneTension_hpp

#include <type_traits>

#include "Mechanics/ForceField/ForceField.h"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "utility.h"

struct MembraneTension : ForceField {
private:

    // Culprit membrane
    const Membrane* membraneCulprit_ = nullptr;

public:

    struct TempMembraneParams {
        struct VertexInfo {
            // the index to vertices in the vectorized coordinate array
            int index = 0;
            bool notBorder = false;
        };
        std::vector< std::array< VertexInfo, 3 >> vertexSet;

        // surface tension
        double tension = 0;
    };

    // (temp) holds vectorized membrane information
    std::vector< TempMembraneParams > tempMembranes;


    virtual void vectorize(const FFCoordinateStartingIndex& si) override {
        using namespace std;
        using AT = Membrane::MeshAttributeType;
        using MT = Membrane::MeshType;

        tempMembranes.clear();

        for(auto m : Membrane::getMembranes()) {

            auto& memParams = tempMembranes.emplace_back();

            const auto& mesh = m->getMesh();
            // Currently applies to general coordinate system
            if(mesh.metaAttribute().vertexSystem == AT::VertexSystem::general) {

                for(const auto& t : mesh.getTriangles()) {
                    const auto vis = medyan::vertexIndices(mesh, t);

                    memParams.vertexSet.push_back({
                        TempMembraneParams::VertexInfo { (int)(mesh.attribute(vis[0]).vertex->getIndex() * 3 + si.vertex), !mesh.isVertexOnBorder(vis[0]) },
                        TempMembraneParams::VertexInfo { (int)(mesh.attribute(vis[1]).vertex->getIndex() * 3 + si.vertex), !mesh.isVertexOnBorder(vis[1]) },
                        TempMembraneParams::VertexInfo { (int)(mesh.attribute(vis[2]).vertex->getIndex() * 3 + si.vertex), !mesh.isVertexOnBorder(vis[2]) }
                    });
                }

                memParams.tension = m->mMembrane.tension;
            }

        }
    }

    virtual floatingpoint computeEnergy(floatingpoint* coord, bool stretched) override {
        using namespace std;

        double en = 0;

        for(const auto& memParams : tempMembranes) {
            double area = 0;

            for(auto& t : memParams.vertexSet) {
                area += medyan::area(
                    makeRefVec<3>(coord + t[0].index),
                    makeRefVec<3>(coord + t[1].index),
                    makeRefVec<3>(coord + t[2].index)
                );
            }

            en += memParams.tension * area;
        }

        return en;
    }
    virtual void computeForces(floatingpoint* coord, floatingpoint* force) override {
        using namespace std;

        for(const auto& memParams : tempMembranes) {
            for(auto& t : memParams.vertexSet) {
                auto rv0 = makeRefVec<3>(coord + t[0].index);
                auto rv1 = makeRefVec<3>(coord + t[1].index);
                auto rv2 = makeRefVec<3>(coord + t[2].index);
                const auto [area, da0, da1, da2] = medyan::areaAndDerivative(rv0, rv1, rv2);

                makeRefVec<3>(force + t[0].index) -= memParams.tension * da0 * t[0].notBorder;
                makeRefVec<3>(force + t[1].index) -= memParams.tension * da1 * t[1].notBorder;
                makeRefVec<3>(force + t[2].index) -= memParams.tension * da2 * t[2].notBorder;
            }
        }
    }

    virtual std::string getName() override { return "Membrane tension"; }

    // Useless overrides
    virtual void whoIsCulprit() override {}
    virtual void cleanup() override {}
    virtual void computeLoadForces() override {}
    virtual std::vector<NeighborList*> getNeighborLists() override { return {}; }
    virtual std::vector<std::string> getinteractionnames() override { return { getName() }; }
};

#endif
