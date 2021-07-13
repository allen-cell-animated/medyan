#ifndef MEDYAN_Mechanics_ForceField_Membrane_MembraneTriangleProtect_hpp
#define MEDYAN_Mechanics_ForceField_Membrane_MembraneTriangleProtect_hpp

#include <array>
#include <type_traits>
#include <vector>

#include "Mechanics/ForceField/ForceField.h"
#include "Mechanics/ForceField/Types.hpp"
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
        const auto deda = k * relDiff / ((1 - relDiff * relDiff) * minArea);

        for(int i = 0; i < forces.size(); ++i) forces[i] -= deda * (*dArea[i]);
    }
};

template< typename Impl, bool warn = false >
class MembraneTriangleProtect: public ForceField {
private:
    using Mems_     = decltype(Membrane::getMembranes());
    using VecDArea_ = decltype(GHalfEdge::dTriangleArea);

    template< typename Float >
    static auto biVec(Float* coordVec, std::array< std::size_t, 3 > bi) {
        using namespace mathfunc;
        using RVT = decltype(makeRefVec<3>(coordVec));

        return std::array< RVT, 3 > {
            makeRefVec< 3 >(coordVec + bi[0]),
            makeRefVec< 3 >(coordVec + bi[1]),
            makeRefVec< 3 >(coordVec + bi[2])
        };
    }

    // Temporary topological information
    //---------------------------------
    Mems_ mems_ = Membrane::getMembranes();

    std::vector< std::array< std::size_t, 3 > >       coordIndices_;
    std::vector< std::array< const VecDArea_*, 3 > >  allDArea_;
    std::vector< const double* >                      allArea_;
    std::vector< const double* >                      allAreaS_;
    std::vector< double >                             initArea_;

public:
    // Parameters
    //---------------------------------
    const double k           = 1000.0;
    const double minRelArea  = 0.35;
    const double warnRelArea = 0.35;

    virtual void vectorize(const FFCoordinateStartingIndex& si) override {
        using MT = Membrane::MeshType;

        coordIndices_.clear();
        allDArea_.clear();
        allArea_.clear();
        allAreaS_.clear();
        initArea_.clear();

        for(const auto m : mems_) {
            const auto& mesh = m->getMesh();
            medyan::assertValidIndexCacheForFF(mesh);

            for(const auto& t : mesh.getTriangles()) {
                coordIndices_.emplace_back();
                allDArea_.emplace_back();
                {
                    auto& biBack = coordIndices_.back();
                    auto& daBack = allDArea_.back();
                    std::size_t bi = 0;
                    mesh.forEachHalfEdgeInPolygon(t, [&](MT::HalfEdgeIndex hei) {
                        biBack[bi] = mesh.attribute(mesh.target(hei)).cachedCoordIndex;
                        daBack[bi] = &mesh.attribute(hei).gHalfEdge.dTriangleArea;
                        ++bi;
                    });
                }
                allArea_.push_back(&t.attr.gTriangle.area);
                allAreaS_.push_back(&t.attr.gTriangleS.area);
                initArea_.push_back(t.attr.gTriangle.area);
            }
        }
    }

    virtual floatingpoint computeEnergy(floatingpoint* coord, bool stretched) override {
        floatingpoint e = 0.0;

        for(std::size_t ti = 0; ti < coordIndices_.size(); ++ti) {

            if(warn && !stretched && *allArea_[ti] < initArea_[ti] * warnRelArea) {
                LOG(WARNING) << "Triangle (" << ti << ") size becomes too small:"
                    << " area=" << *allArea_[ti] << " init_area=" << initArea_[ti];
            }

            e += Impl::energy(*(stretched ? allAreaS_ : allArea_)[ti] / (initArea_[ti] * minRelArea), k);
        }

        return e;
    }
    virtual void computeForces(floatingpoint* coord, floatingpoint* force) override {
        for(std::size_t ti = 0; ti < coordIndices_.size(); ++ti) {

            if(warn && *allArea_[ti] < initArea_[ti] * warnRelArea) {
                LOG(WARNING) << "Triangle (" << ti << ") size becomes too small:"
                    << " area=" << *allArea_[ti] << " init_area=" << initArea_[ti];
            }

            Impl::force(
                biVec(force, coordIndices_[ti]),
                *allArea_[ti], initArea_[ti] * minRelArea, allDArea_[ti], k
            );
        }
    }

    virtual string getName() override { return "Membrane Triangle Protect"; }

    // Useless overrides
    virtual void cleanup() override {}
    virtual void computeLoadForces() override {}
    virtual void whoIsCulprit() override {}
    virtual std::vector<NeighborList*> getNeighborLists() override { return std::vector<NeighborList*>(); }
    virtual std::vector<std::string> getinteractionnames() override { return { getName() }; }
};


#endif
