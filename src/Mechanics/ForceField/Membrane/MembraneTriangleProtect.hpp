#ifndef MEDYAN_Mechanics_ForceField_Membrane_MembraneTriangleProtect_hpp
#define MEDYAN_Mechanics_ForceField_Membrane_MembraneTriangleProtect_hpp

#include <array>
#include <type_traits>
#include <vector>

#include "Mechanics/ForceField/Membrane/MembraneInteractions.hpp"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Util/Io/Log.hpp"
#include "Util/Math/Vec.hpp"

struct MembraneTriangleProtectFene {

    static floatingpoint energy(double relArea, double k) {
        if(relArea >= 1.0) return 0.0;

        const auto relDiff = relArea - 1;
        return -0.5 * k * std::log(1 - relDiff * relDiff);
    }

    template< typename RefVecType, typename VecType >
    static void force(
        std::array< RefVecType, 3 > forces, 
        double area, double minArea, const std::array< VecType*, 3 >& dArea, double k
    ) {
        if(area >= minArea) return;

        const auto relDiff = area / minArea - 1;
        const auto deda = -k * relDiff / ((1 - relDiff * relDiff) * minArea);

        for(int i = 0; i < forces.size(); ++i) forces[i] -= deda * (*dArea[i]);
    }
};

template< typename Impl, bool warn = false >
class MembraneTriangleProtect: public MembraneInteractions {
private:
    using Mems_     = decltype(Membrane::getMembranes());
    using VecDArea_ = decltype(GHalfEdge::dTriangleArea);

    template< typename Float >
    static auto biVec(Float* coordVec, std::array< std::size_t, 3 > bi) {
        using namespace mathfunc;
        using RVT = decltype(makeRefVec<3>(coordVec));

        return std::array< RVT, 3 > {
            makeRefVec< 3 >(coordVec + 3 * bi[0]),
            makeRefVec< 3 >(coordVec + 3 * bi[1]),
            makeRefVec< 3 >(coordVec + 3 * bi[2])
        };
    }

    // Temporary topological information
    //---------------------------------
    Mems_ mems_ = Membrane::getMembranes();

    std::vector< std::array< std::size_t, 3 > >       beadIndices_;
    std::vector< std::array< const VecDArea_*, 3 > >  allDArea_;
    std::vector< const double* >                      allArea_;
    std::vector< const double* >                      allAreaS_;
    std::vector< double >                             initArea_;

public:
    // Parameters
    //---------------------------------
    const double k           = 1.0;
    const double minRelArea  = 0.1;
    const double warnRelArea = 0.05;

    virtual void vectorize() override {
        beadIndices_.clear();
        allDArea_.clear();
        allArea_.clear();
        allAreaS_.clear();
        initArea_.clear();

        for(const auto m : mems_) {
            const auto& mesh = m->getMesh();
            for(const auto& t : mesh.getTriangles()) {
                beadIndices_.emplace_back();
                allDArea_.emplace_back();
                {
                    auto& biBack = beadIndices_.back();
                    auto& daBack = allDArea_.back();
                    std::size_t bi = 0;
                    mesh.forEachHalfEdgeInPolygon(t, [&](std::size_t hei) {
                        biBack[bi] = mesh.getVertexAttribute(mesh.target(hei)).vertex->getStableIndex();
                        daBack[bi] = &mesh.getHalfEdgeAttribute(hei).gHalfEdge.dTriangleArea;
                        ++bi;
                    });
                }
                allArea_.push_back(&t.attr.gTriangle.area);
                allAreaS_.push_back(&t.attr.gTriangleS.area);
                initArea_.push_back(t.attr.gTriangle.area);
            }
        }
    }

    virtual floatingpoint computeEnergy(const floatingpoint* coord, bool stretched) override {
        floatingpoint e = 0.0;

        for(std::size_t ti = 0; ti < beadIndices_.size(); ++ti) {

            if(warn && !stretched && *allArea_[ti] < initArea_[ti] * warnRelArea) {
                LOG(WARNING) << "Triangle (" << ti << ") size becomes too small:"
                    << " area=" << *allArea_[ti] << " init_area=" << initArea_[ti];
            }

            e += Impl::energy(*(stretched ? allAreaS_ : allArea_)[ti] / (initArea_[ti] * minRelArea), k);
        }

        return e;
    }
    virtual void computeForces(const floatingpoint* coord, floatingpoint* force) override {
        for(std::size_t ti = 0; ti < beadIndices_.size(); ++ti) {

            if(warn && *allArea_[ti] < initArea_[ti] * warnRelArea) {
                LOG(WARNING) << "Triangle (" << ti << ") size becomes too small:"
                    << " area=" << *allArea_[ti] << " init_area=" << initArea_[ti];
            }

            Impl::force(
                biVec(force, beadIndices_[ti]),
                *allArea_[ti], initArea_[ti] * minRelArea, allDArea_[ti], k
            );
        }
    }

    virtual string getName() const override { return "Membrane Triangle Protect"; }

};


#endif
