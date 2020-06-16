// WARNING: THIS FILE IS NOT USABLE FOR NOW

#ifndef MEDYAN_Mechanics_ForceField_Membrane_MembraneTension_hpp
#define MEDYAN_Mechanics_ForceField_Membrane_MembraneTension_hpp

#include <type_traits>

#include "Mechanics/ForceField/ForceField.h"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "utility.h"

class MembraneTension: public ForceField {
private:

    // Culprit membrane
    const Membrane* membraneCulprit_ = nullptr;

public:
    virtual floatingpoint computeEnergy(floatingpoint* coord, bool stretched) override {
        using namespace std;

        double U = 0;

        for(auto m: Membrane::getMembranes()) {

            double enMem = 0.0;

            if constexpr(type == MembraneStretchingType::localHarmonic) {
            }
            else {
                double area = 0.0;

                for(const auto& t : m->getMesh().getTriangles())
                    area += stretched ? t.attr.gTriangleS.area : t.attr.gTriangle.area;

                enMem = MembraneStretchingHarmonic{}.energy(
                    area,
                    m->getMMembrane()->getKElastic();
                    m->getMMembrane()->getEqArea();
                );

            IF_CONSTEXPR (std::is_same< Impl, MembraneStretchingHarmonic >::value) {
                const auto kElastic = m->getMMembrane()->getKElastic();
                const auto eqArea = m->getMMembrane()->getEqArea();

    #ifdef __cpp_if_constexpr
                U_i = _impl.energy(area, kElastic, eqArea);
    #else
                U_i = implEnergy(_impl, area, kElastic, eqArea);
    #endif

            } else {
                const auto tension = m->getMMembrane()->getTension();

    #ifdef __cpp_if_constexpr
                U_i = _impl.energy(area, tension);
    #else
                U_i = implEnergy(_impl, area, tension, 0.0);
    #endif
            }

            if(fabs(U_i) == numeric_limits<double>::infinity()
            || U_i != U_i || U_i < -1.0) { // (U_i != U_i) => (U_i == NaN)
                
                //set culprit and return
                membraneCulprit_ = m;
                
                return -1;
            }
            else
                U += U_i;

        }

        return U;

    }
    virtual void computeForces(floatingpoint* coord, floatingpoint* force) override;

    virtual string getName() override { return "Membrane tension"; }

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
    virtual std::vector<NeighborList*> getNeighborLists() override { return std::vector<NeighborList*>(); }
};

#endif
