#ifndef MEDYAN_Mechanics_ForceField_Membrane_MembraneStretching_hpp
#define MEDYAN_Mechanics_ForceField_Membrane_MembraneStretching_hpp

#include <limits>
#include <type_traits>

#include "Mechanics/ForceField/ForceField.h"
#include "Mechanics/ForceField/Membrane/MembraneStretchingImpl.hpp"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "utility.h"

enum class MembraneStretchingType {
    localHarmonic,
    globalHarmonic
};

template< MembraneStretchingType type >
class MembraneStretching: public ForceField {
private:
    // Culprit membrane
    const Membrane* membraneCulprit_ = nullptr;

public:
    virtual floatingpoint computeEnergy(floatingpoint* coord, bool stretched) override {
        using namespace std;

        double en = 0;

        for(auto m: Membrane::getMembranes()) {

            double enMem = 0.0;

            if constexpr(type == MembraneStretchingType::localHarmonic) {
                for(const auto& t : m->getMesh().getTriangles()) {
                    const auto& area = stretched ? t.attr.gTriangleS.area : t.attr.gTriangle.area;
                    enMem += MembraneStretchingHarmonic{}.energy(
                        area,
                        t.attr.triangle->mTriangle.kArea,
                        t.attr.triangle->mTriangle.eqArea
                    )
                }
            }
            else {
                double area = 0.0;

                for(const auto& t : m->getMesh().getTriangles())
                    area += stretched ? t.attr.gTriangleS.area : t.attr.gTriangle.area;

                enMem = MembraneStretchingHarmonic{}.energy(
                    area,
                    m->getMMembrane()->getKElastic(),
                    m->getMMembrane()->getEqArea()
                );
            }

            if(!isfinite(enMem)) { // (U_i != U_i) => (U_i == NaN)
                
                //set culprit and return
                membraneCulprit_ = m;
                
                return numeric_limits<double>::infinity();
            }
            else
                en += enMem;

        }

        return en;
    }

    virtual void computeForces(floatingpoint* coord, floatingpoint* force) override {
        using namespace std;

        // Configure force buffer
        constexpr bool useForceBuffer = false;
        const std::size_t dof = 0; // TODO get dof

        if (useForceBuffer) {
            forceBuffer_.assign(dof, 0.0);
        }
        floatingpoint* const f = useForceBuffer ? forceBuffer_.data() : force;

        for (auto m: Membrane::getMembranes()) {

            Membrane::MembraneMeshAttributeType::cacheIndices(m->getMesh());

            const auto& mesh = m->getMesh();

            if constexpr(type == MembraneStretchingType::localHarmonic) {
                const size_t numTriangles = mesh.getTriangles().size();
                for(size_t ti = 0; ti < numTriangles; ++ti) {
                    const auto& ta = mesh.getTriangleAttribute(ti);

                    for(size_t i = 0; i < 3; ++i) {
                        const size_t hei = ta.cachedHalfEdgeIndex[i];
                        const auto& dArea = mesh.getHalfEdgeAttribute(hei).gHalfEdge.dTriangleArea;

                        MembraneStretchingHarmonic{}.forces(
                            f + 3 * ta.cachedCoordIndex[i],
                            ta.gTriangle.area, dArea,
                            ta.triangle->mTriangle.kArea, ta.triangle->mTriangle.eqArea
                        );
                    }
                }

            }
            else {

                const auto kElastic = m->getMMembrane()->getKElastic();
                const auto eqArea = m->getMMembrane()->getEqArea();

                double area = 0.0;
                for(const auto& t : mesh.getTriangles()) area += t.attr.gTriangle.area;

                const size_t numTriangles = mesh.getTriangles().size();
                for(size_t ti = 0; ti < numTriangles; ++ti) {
                    const auto& ta = mesh.getTriangleAttribute(ti);

                    for(size_t i = 0; i < 3; ++i) {
                        const size_t hei = ta.cachedHalfEdgeIndex[i];
                        const auto& dArea = mesh.getHalfEdgeAttribute(hei).gHalfEdge.dTriangleArea;

                        MembraneStretchingHarmonic{}.forces(
                            f + 3 * ta.cachedCoordIndex[i],
                            area, dArea, kElastic, eqArea
                        );
                    }
                }
            } // end if (type == ...)

        } // end for (m : membranes)

        if(useForceBuffer) {
            std::transform(
                force, force + dof, forceBuffer_.begin(),
                force,
                std::plus<>{}
            );
        }

    }

    virtual string getName() override { return "Membrane Stretching"; }

    virtual void whoIsCulprit() override {
        if(membraneCulprit_) {
            LOG(INFO) << "Printing culprit membrane...";
            membraneCulprit_->printSelf();
        } else {
            LOG(ERROR) << "Membrane culprit is not set.";
        }
    }

    // Useless overrides
    virtual void vectorize() override {}
    virtual void cleanup() override {}
    virtual void computeLoadForces() override {}
    virtual std::vector<NeighborList*> getNeighborLists() override { return {}; }
};

#endif
