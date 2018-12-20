#ifndef MEDYAN_Membrane_h
#define MEDYAN_Membrane_h

#include <array>
#include <tuple>
#include <vector>
#include <memory>

#include "Database.h"
#include "Geometric.h"
#include "Trackable.h"
#include "Composite.h"

#include "GMembrane.h"
#include "MMembrane.h"
#include "Structure/SurfaceMesh/Edge.h"
#include "Structure/SurfaceMesh/SurfaceMesh.hpp"
#include "Structure/SurfaceMesh/Triangle.h"
#include "Structure/SurfaceMesh/Vertex.h"

// FORWARD DECLARATIONS
class SubSystem;

// Element data specification
struct HalfEdgeData {};
using MembraneMeshData = SurfaceTriangularMesh::TriangularMeshAttribute< VertexData, EdgeData, HalfEdgeData, TriangleData >;

/******************************************************************************
Topologically, the membrane is represented by a 2d surface with 2 sides (which
means no Klein bottles are allowed!). The surface is constructed by
interconnected vertices, edges and triangles.

The Membrane class is a container holding all the relative vertices, edges and
triangles. Meshwork initialization and geometry update are all managed by this
class.
******************************************************************************/

class Membrane: public Composite, public Trackable, public Geometric {

private:

    SurfaceTriangularMeshData<MembraneMeshData> _mesh;

    vector<Triangle*> _triangleVector; // collection of triangles
    vector<Edge*> _edgeVector; // collection of edges
    vector<Vertex*> _vertexVector; // collection of vertices

    unique_ptr<GMembrane> _gMembrane; // pointer to geometric membrane object
    unique_ptr<MMembrane> _mMembrane; // pointer to mechanical membrane object

    short _memType; // Membrane type

    SubSystem* _subSystem; // SubSystem pointer

    static Database<Membrane*> _membranes; // Collection in SubSystem
    int _id; // Unique integer id of this membrane

    bool _isClosed = true; // Whether the membrane is topologically closed, regardless of genus
    int _genus = 0; // Genus of the surface. Normally 0, as for a topologically spherical shape

public:
    // Constructors
    // This constructor creates a membrane according to vertex and neighbor data
    Membrane(SubSystem* s, short membraneType,
        const vector<tuple<array<double, 3>, vector<size_t>>>& membraneData);

    // This destructor is called when a membrane is to be removed from the system.
    // Removes all vertices, edges and triangles associated with the membrane.
    ~Membrane();

    /// Get vector of triangles/edges/vertices that this membrane contains.
    vector<Triangle*>& getTriangleVector() {return _triangleVector;}
    vector<Edge*>& getEdgeVector() { return _edgeVector; }
    vector<Vertex*>& getVertexVector() { return _vertexVector; }

    // Get Id
    int getId()const { return _id; }
    
    
    // SubSystem management, inherited from Trackable
    virtual void addToSubSystem()override { _membranes.addElement(this); }
    virtual void removeFromSubSystem()override { _membranes.removeElement(this); }
    
    /// Get all instances of this class from the SubSystem
    static const vector<Membrane*>& getMembranes() {
        return _membranes.getElements();
    }
    /// Get the number of membranes in this system
    static int numMembranes() {
        return _membranes.countElements();
    }

    //@{
    /// Implements Component
    // Get type
    int getType()override { return _memType; }
    // Print self information
    virtual void printSelf()const override;
    //@}

    /**************************************************************************
    Geometric
    **************************************************************************/
    /**
     * Implements Geometric.
     * If calcDerivative is true, most implementation would assume d is zero
     *   regardless of the actual passed d value. If calcDerivative is false,
     *   most implementation would store the result in member variables with
     *   "stretched" in their name.
     */
    virtual void updateGeometry(bool calcDerivative=false, double d=0.0)override;

    // Get geo membrane
    GMembrane* getGMembrane() { return _gMembrane.get(); }

    /**
     * Use pseudo normal signed distance field method to get the signed distance to a point.
     * If the point is outside, the result is positive and vice versa.
     * Throws an exception if the membrane is not closed.
     * The function will search through the whole meshwork, so it might not be efficient.
     * However, if the "safe" flag is turned off and neighboring compartments happen to contain membrane elements,
     *   search space will be limited to those compartments to save time. In this case meshwork size should be
     *   much smaller than the compartment size.
     */
    double signedDistance(const std::array<double, 3>& p, bool safe=false)const;
    /**
     * Use signed distance or other methods to judge whether a point is inside membrane.
     * Throws an exception if the membrane is not closed.
     */
    bool contains(const std::array<double, 3>& point)const;

    // Function to monitor the quality of the meshwork
    double meshworkQuality()const; // Must be used after updating the geometry
                                   // Returns a value between 0 and 1,
                                   // 1 being best and 0 being worst.

    /**************************************************************************
    Topological
    **************************************************************************/
    bool isClosed()const { return _isClosed; }

    /**************************************************************************
    Mechanics
    **************************************************************************/
    // Get mech membrane
    MMembrane* getMMembrane() { return _mMembrane.get(); }


};





#endif
