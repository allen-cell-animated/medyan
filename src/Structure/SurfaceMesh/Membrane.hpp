#ifndef MEDYAN_Membrane_hpp
#define MEDYAN_Membrane_hpp

#include <array>
#include <limits> // numeric_limits
#include <memory>
#include <stdexcept> // runtime_error
#include <vector>

#include "Composite.h"
#include "MathFunctions.h"
#include "Structure/Database.h"
#include "Structure/SubSystem.h"
#include "Structure/SurfaceMesh/MembraneHierarchy.hpp"
#include "Structure/SurfaceMesh/MembraneMeshAttribute.hpp"
#include "Structure/SurfaceMesh/MMembrane.hpp"
#include "Structure/SurfaceMesh/SurfaceMesh.hpp"
#include "Structure/Trackable.h"
#include "SysParams.h"

/******************************************************************************
Topologically, the membrane is represented by a 2d surface with 2 sides (which
means no Klein bottles are allowed!). The surface is constructed by
interconnected vertices, edges and triangles.

The Membrane class is a manager for constructing the mesh and computing the
geometry. It contains a meshwork instance that's responsible for adding and
removing vertices, edges (halfedges) and triangles to/from the SubSystem.
However, the ownership of all elements is in this Membrane class through
inheriting Composite.
******************************************************************************/
class Membrane: public Composite, public Trackable, public Database< Membrane, false > {
public:
    using MembraneMeshAttributeType = MembraneMeshAttribute< HalfEdgeMesh >;
    using coordinate_type = typename MembraneMeshAttributeType::coordinate_type;
    using MeshType = HalfEdgeMesh< MembraneMeshAttributeType >;

    using HierarchyType = MembraneHierarchy< Membrane >;

private:

    MeshType mesh_;

    unique_ptr<MMembrane> mMembrane_; // pointer to mechanical membrane object

    short memType_; // Membrane type

    SubSystem* _subSystem; // SubSystem pointer

public:

    // Constructors
    // This constructor creates a membrane according to vertex and neighbor data
    Membrane(
        SubSystem* s,
        short membraneType,
        const std::vector< coordinate_type >& vertexCoordinateList,
        const std::vector< std::array< size_t, 3 > >& triangleVertexIndexList
    ) : Trackable(false, false, false, false),
        mesh_(typename MembraneMeshAttributeType::MetaAttribute{s, this}),
        _subSystem(s), memType_(membraneType) {
        
        // Build the meshwork topology using vertex and triangle information
        mesh_.init<typename MeshType::VertexTriangleInitializer>(
            vertexCoordinateList.size(),
            triangleVertexIndexList,
            typename MembraneMeshAttributeType::AttributeInitializerInfo{ vertexCoordinateList }
        );

        // Update geometry
        updateGeometryValue<false>();

        initMMembrane();

        // Add to membrane hierarchy (must have updated geometry)
        HierarchyType::addMembrane(this);
    }

    ~Membrane() {
        HierarchyType::removeMembrane(this);
    }

    /// Get vector of triangles/edges/vertices that this membrane contains.
    const auto& getMesh() const { return mesh_; }
    auto&       getMesh()       { return mesh_; }

    // Helper function to initialize MMembrane
    void initMMembrane() {
        // Calculate the total area and volume to set the equilibrium area and volume
        double area = 0.0;
        double volume = 0.0;
        for(const auto& t : mesh_.getTriangles()) {
            area += t.attr.gTriangle.area;
            volume += t.attr.gTriangle.coneVolume;
        }

        mMembrane_ = std::make_unique<MMembrane>(memType_);
        mMembrane_->setEqArea(
            area * SysParams::Mechanics().MemEqAreaFactor[memType_]
        );
        mMembrane_->setEqVolume(volume);
    } // void initMMembrane()

    // SubSystem management, inherited from Trackable
    virtual void addToSubSystem()override { }
    virtual void removeFromSubSystem()override { }
    
    /// Get all instances of this class from the SubSystem
    static const vector<Membrane*>& getMembranes() {
        return getElements();
    }
    /// Get the number of membranes in this system
    static std::size_t numMembranes() {
        return getElements().size();
    }

    /// Implements Component
    // Get type
    int getType() override { return memType_; }
    // Print self information
    virtual void printSelf() const override {
        using namespace std;

        cout << endl;

        cout << "Membrane Id = " << getId() << endl;
        cout << "Membrane type = " << memType_ << endl;
        cout << "Number of vertices, edges, half edges, triangles, borders =\n  "
            << mesh_.numVertices() << ' ' << mesh_.numEdges() << ' ' << mesh_.numHalfEdges() << ' '
            << mesh_.numTriangles() << ' ' << mesh_.numBorders() << endl;

        cout << endl;
    }

    /**************************************************************************
    Geometric
    **************************************************************************/
    template<
        bool stretched = false,
        MembraneMeshAttributeType::GeometryCurvaturePolicy curvPol
            = MembraneMeshAttributeType::GeometryCurvaturePolicy::withSign
    > void updateGeometryValue() {
        MembraneMeshAttributeType::template updateGeometryValue<stretched, curvPol>(mesh_);
    }
    template<
        MembraneMeshAttributeType::GeometryCurvaturePolicy curvPol
            = MembraneMeshAttributeType::GeometryCurvaturePolicy::withSign
    > void updateGeometryValueWithDerivative() {
        MembraneMeshAttributeType::template updateGeometryValueWithDerivative<curvPol>(mesh_);
    }
    void updateGeometryValueForSystem() {
        MembraneMeshAttributeType::updateGeometryValueForSystem(mesh_);
    }

    /**
     * Use pseudo normal signed distance field method to get the signed distance to a point.
     * If the point is outside, the result is positive and vice versa.
     * Throws an exception if the membrane is not closed.
     * The function will search through the whole meshwork, so it might not be efficient.
     */
    template< typename VecType, std::enable_if_t< VecType::vec_size == 3 >* = nullptr >
    double signedDistance(const VecType& p, bool allowOpen = false) const {
        if(!allowOpen && !isClosed())
            throw std::runtime_error("Membrane is not closed while trying to find signed distance field.");
        return MembraneMeshAttributeType::signedDistance(mesh_, p);
    }
    /**
     * Use signed distance or other methods to judge whether a point is inside membrane.
     * Throws an exception if the membrane is not closed.
     */
    template< typename VecType, std::enable_if_t< VecType::vec_size == 3 >* = nullptr >
    bool contains(const VecType& p) const {
        if(!isClosed()) throw std::runtime_error("Membrane is not closed while trying to find signed distance field.");
        return MembraneMeshAttributeType::contains(mesh_, p);
    }

    /**************************************************************************
    Topological
    **************************************************************************/
    bool isClosed() const { return mesh_.isClosed(); }

    /**************************************************************************
    Mechanics
    **************************************************************************/
    // Get mech membrane
    MMembrane* getMMembrane() { return mMembrane_.get(); }


};



#endif
