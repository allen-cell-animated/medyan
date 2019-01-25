#ifndef MEDYAN_SurfaceMesh_hpp
#define MEDYAN_SurfaceMesh_hpp

#include <algorithm>
#include <array>
#include <cstddef> // ptrdiff_t
#include <iterator>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

/******************************************************************************
The data structure for an orientable, manifold 2d triangular meshwork in 3d
space. The data structure only provides topological relationship, and is
completed by the Attribute class.

The connectivity is similar to CGAL library, which is halfedge based.

target:     The vertex that the halfedge points to.
opposite:   The halfedge on the same edge but in the opposite direction.
face:       The triangle that the halfedge is circling ccw.
next:       The next halfedge in the current face circling ccw.
prev:       The previous halfedge in the current face circling ccw.
edge:       The undirected edge of this halfedge.

All the other elements must have at least one index pointing to a halfedge.
******************************************************************************/

// The DeletableVector is a thin wrapper for std::vector.
// This container is designed for trivially_copyable types.
// When an element is removed, instead of doing vector::erase,
// it essentially swaps the element with the last one, and pops the vector.
template< typename T > class DeletableVector {
public:
    struct IndexMove {
        size_t from, to;
        bool valid;
    };

private:
    std::vector< T > _value;

public:

    using iterator       = typename std::vector< T >::iterator;
    using const_iterator = typename std::vector< T >::const_iterator;
    iterator       begin() noexcept       { return _value.begin(); }
    const_iterator begin() const noexcept { return _value.begin(); }
    iterator       end() noexcept       { return _value.end(); }
    const_iterator end() const noexcept { return _value.end(); }

    // Insert a new element. Returns the new index.
    size_t insert() {
        _value.emplace_back();
        return _value.size() - 1;
    }

    // Remove an element from the container.
    // Might change the position of certain elements, indicated by the return value.
    // If index is out of range, the behavior is undefined.
    IndexMove erase(size_t index) {
        IndexMove res;
        const size_t lastIndex = _value.size() - 1;
        if(index == lastIndex) {
            _value.pop_back();
            res.valid = false;
        } else {
            // Move value from lastIndex to index
            _value[index] = _value[lastIndex];
            _value.pop_back();
            res.from = lastIndex;
            res.to = index;
            res.valid = true;
        }
        return res;
    }

    size_t size() const noexcept { return _value.size(); }

    T&       operator[](size_t index)       { return _value[index]; }
    const T& operator[](size_t index) const { return _value[index]; }

    // This function should only be called during initialization/finalization
    std::vector< T >& getValue() { return _value; }

};

// The Attribute class must implement
//   - Type VertexAttribute
//     - void setIndex(size_t)
//   - Type EdgeAttribute
//     - void setIndex(size_t)
//   - Type HalfEdgeAttribute
//     - void setIndex(size_t)
//   - Type TriangleAttribute
//     - void setIndex(size_t)
//   - Type MetaAttribute
//   - Type AttributeInitializerInfo
//   - void init(Mesh, const AttributeInitializerInfo&)
//   - AttributeInitializerInfo extract(Mesh)
template< typename Attribute > class SurfaceTriangularMesh {
public:

    using AttributeType = Attribute;
    using MeshType = SurfaceTriangularMesh;

    using VertexAttribute   = typename Attribute::VertexAttribute;
    using EdgeAttribute     = typename Attribute::EdgeAttribute;
    using HalfEdgeAttribute = typename Attribute::HalfEdgeAttribute;
    using TriangleAttribute = typename Attribute::TriangleAttribute;
    using MetaAttribute     = typename Attribute::MetaAttribute;

    // The elements should be trivially copyable.
    struct Vertex {
        size_t halfEdgeIndex; // Only one HalfEdge targeting the vertex is needed.
        size_t degree; // Number of neighbors
        VertexAttribute attr;
    };
    struct HalfEdge {
        bool hasOpposite;
        size_t triangleIndex;
        size_t targetVertexIndex;
        size_t oppositeHalfEdgeIndex;
        size_t nextHalfEdgeIndex;
        size_t prevHalfEdgeIndex;
        size_t edgeIndex;
        HalfEdgeAttribute attr;
    };
    struct Edge {
        size_t halfEdgeIndex; // Only one HalfEdge is needed.
        EdgeAttribute attr;
    };
    struct Triangle {
        size_t halfEdgeIndex; // Only one HalfEdge is needed.
        TriangleAttribute attr;
    };

private:

    DeletableVector<Triangle> _triangles; // collection of triangles
    DeletableVector<HalfEdge> _halfEdges; // collection of halfedges
    DeletableVector<Edge>     _edges;     // collection of edges
    DeletableVector<Vertex>   _vertices;  // collection of vertices

    MetaAttribute _meta;

    bool _isClosed;
    int _genus = 0; // Genus of the surface. Currently it is not tracked.

    // Element accessor
    template< typename Element, std::enable_if_t<std::is_same<Element, Triangle>::value, void>* = nullptr>
    auto& _getElements() { return _triangles; }
    template< typename Element, std::enable_if_t<std::is_same<Element, HalfEdge>::value, void>* = nullptr>
    auto& _getElements() { return _halfEdges; }
    template< typename Element, std::enable_if_t<std::is_same<Element, Edge>::value, void>* = nullptr>
    auto& _getElements() { return _edges; }
    template< typename Element, std::enable_if_t<std::is_same<Element, Vertex>::value, void>* = nullptr>
    auto& _getElements() { return _vertices; }

    // Meshwork registration helper
    void _registerTriangle(size_t ti, size_t hei0, size_t hei1, size_t hei2) {
        _triangles[ti].halfEdgeIndex = hei0;
        _halfEdges[hei0].nextHalfEdgeIndex = hei1;
        _halfEdges[hei0].prevHalfEdgeIndex = hei2;
        _halfEdges[hei0].triangleIndex = ti;
        _halfEdges[hei1].nextHalfEdgeIndex = hei2;
        _halfEdges[hei1].prevHalfEdgeIndex = hei0;
        _halfEdges[hei1].triangleIndex = ti;
        _halfEdges[hei2].nextHalfEdgeIndex = hei0;
        _halfEdges[hei2].prevHalfEdgeIndex = hei1;
        _halfEdges[hei2].triangleIndex = ti;
    }
    void _registerEdge(size_t ei, size_t hei0, size_t hei1) {
        _edges[ei].halfEdgeIndex = hei0;
        _halfEdges[hei0].hasOpposite = true;
        _halfEdges[hei0].oppositeHalfEdgeIndex = hei1;
        _halfEdges[hei0].edgeIndex = ei;
        _halfEdges[hei1].hasOpposite = true;
        _halfEdges[hei1].oppositeHalfEdgeIndex = hei0;
        _halfEdges[hei1].edgeIndex = ei;
    }

    template< typename Operation > size_t _newVertex(const Operation& op) {
        size_t index = _vertices.insert();
        Attribute::newVertex(*this, index, op);
        return index;
    }
    template< typename Operation > size_t _newEdge(const Operation& op) {
        size_t index = _edges.insert();
        Attribute::newEdge(*this, index, op);
        return index;
    }
    template< typename Operation > size_t _newHalfEdge(const Operation& op) {
        size_t index = _halfEdges.insert();
        Attribute::newHalfEdge(*this, index, op);
        return index;
    }
    template< typename Operation > size_t _newTriangle(const Operation& op) {
        size_t index = _triangles.insert();
        Attribute::newTriangle(*this, index, op);
        return index;
    }

    void _removeVertex(size_t index) {
        Attribute::removeVertex(*this, index);
        auto moveIndex = _vertices.erase(index);
        if(moveIndex.valid) {
            // Need to update all stored indices/reference/pointer to the vertex.
            forEachHalfEdgeTargetingVertex(index, [this, index](size_t hei) {
                _halfEdges[hei].targetVertexIndex = index;
            });
            _vertices[index].attr.setIndex(index);
        }
    }
    void _removeHalfEdge(size_t index) {
        Attribute::removeHalfEdge(*this, index);
        auto moveIndex = _halfEdges.erase(index);
        if(moveIndex.valid) {
            if(hasOpposite(index)) _halfEdges[opposite(index)].oppositeHalfEdgeIndex = index;
            if(_triangles[triangle(index)].halfEdgeIndex == moveIndex.from)
                _triangles[triangle(index)].halfEdgeIndex = index;
            if(_vertices[target(index)].halfEdgeIndex == moveIndex.from)
                _vertices[target(index)].halfEdgeIndex = index;
            if(_edges[edge(index)].halfEdgeIndex == moveIndex.from)
                _edges[edge(index)].halfEdgeIndex = index;
            _halfEdges[next(index)].prevHalfEdgeIndex = index;
            _halfEdges[prev(index)].nextHalfEdgeIndex = index;
            _halfEdges[index].attr.setIndex(index);
        }
    }
    void _removeEdge(size_t index) {
        Attribute::removeEdge(*this, index);
        auto moveIndex = _edges.erase(index);
        if(moveIndex.valid) {
            forEachHalfEdgeInEdge(index, [this, index](size_t hei) {
                _halfEdges[hei].edgeIndex = index;
            });
            _edges[index].attr.setIndex(index);
        }
    }
    void _removeTriangle(size_t index) {
        Attribute::removeTriangle(*this, index);
        auto moveIndex = _triangles.erase(index);
        if(moveIndex.valid) {
            forEachHalfEdgeInTriangle(index, [this, index](size_t hei) {
                _halfEdges[hei].triangleIndex = index;
            });
            _triangles[index].attr.setIndex(index);
        }
    }

    template< typename Element > void _clearElement() {
        auto& elements = _getElements< Element >();
        for(size_t i = 0; i < elements.size(); ++i)
            Attribute::template removeElement<MeshType, Element>(*this, i);
        elements.getValue().clear();
    }

    void _clear() {
        _clearElement<Vertex>();
        _clearElement<HalfEdge>();
        _clearElement<Edge>();
        _clearElement<Triangle>();
    }

public:

    // Constructors
    SurfaceTriangularMesh(const MetaAttribute& meta) : _meta(meta) {}

    // Destructor
    ~SurfaceTriangularMesh() {
        _clear();
    }

    struct VertexTriangleInitializer {
        struct Info {
            size_t numVertices;
            std::vector< std::array< size_t, 3 > > triangleVertexIndexList;
            typename Attribute::AttributeInitializerInfo attributeInitializerInfo;
        };

        void init(
            SurfaceTriangularMesh& mesh,
            size_t numVertices,
            const std::vector< std::array< size_t, 3 > >& triangleVertexIndexList,
            const typename Attribute::AttributeInitializerInfo& attributeInitializerInfo
        ) const {
            mesh._vertices.getValue().resize(numVertices);
            const size_t numTriangles = triangleVertexIndexList.size();
            mesh._triangles.getValue().resize(numTriangles);
            const size_t numHalfEdges = 3 * numTriangles;
            mesh._halfEdges.getValue().reserve(numHalfEdges);
            mesh._edges.getValue().reserve(numHalfEdges / 2); // Might be more than this number with borders.

            struct VertexAdditionalInfo {
                bool hasTargetingHalfEdge = false;
                std::vector< size_t > leavingHalfEdgeIndices;
            };
            std::vector< VertexAdditionalInfo > vai(numVertices);

            for(size_t ti = 0; ti < numTriangles; ++ti) {
                const auto& t = triangleVertexIndexList[ti];
                mesh._triangles[ti].halfEdgeIndex = mesh._halfEdges.size(); // The next inserted halfedge index
                for(size_t i = 0; i < 3; ++i) {
                    const size_t hei = mesh._halfEdges.insert();
                    HalfEdge& he = mesh._halfEdges[hei];
                    he.hasOpposite = false;
                    he.triangleIndex = ti;
                    he.targetVertexIndex = t[i];
                    he.nextHalfEdgeIndex = (i == 2 ? hei - 2 : hei + 1);
                    he.prevHalfEdgeIndex = (i == 0 ? hei + 2 : hei - 1);

                    // Remembering this edge in the vertices.
                    const size_t leftVertexIndex = t[i == 0 ? 2 : i - 1];
                    vai[leftVertexIndex].leavingHalfEdgeIndices.push_back(hei);
                    // Search in the target vertex, whether there's an opposite halfedge leaving
                    const auto findRes = std::find_if(
                        vai[t[i]].leavingHalfEdgeIndices.begin(),
                        vai[t[i]].leavingHalfEdgeIndices.end(),
                        [&mesh, leftVertexIndex](size_t leavingHalfEdgeIndex) {
                            return leftVertexIndex == mesh._halfEdges[leavingHalfEdgeIndex].targetVertexIndex;
                        }
                    );
                    if(findRes == vai[t[i]].leavingHalfEdgeIndices.end()) {
                        // opposite not found
                        mesh._edges[mesh._edges.insert()].halfEdgeIndex = hei;
                        he.edgeIndex = mesh._edges.size() - 1;
                    } else {
                        // opposite found
                        he.hasOpposite = true;
                        he.oppositeHalfEdgeIndex = *findRes;
                        he.edgeIndex = mesh._halfEdges[he.oppositeHalfEdgeIndex].edgeIndex;

                        mesh._halfEdges[he.oppositeHalfEdgeIndex].hasOpposite = true;
                        mesh._halfEdges[he.oppositeHalfEdgeIndex].oppositeHalfEdgeIndex = hei;
                    }

                    // Set vertex half edge index
                    if(!vai[t[i]].hasTargetingHalfEdge) {
                        mesh._vertices[t[i]].halfEdgeIndex = hei;
                        vai[t[i]].hasTargetingHalfEdge = true;
                    }
                } // end loop halfedges
            } // end loop triangles

            // Registering vertex degrees (not accurate when there are holes)
            for(size_t vi = 0; vi < numVertices; ++vi) {
                mesh._vertices[vi].degree = vai[vi].leavingHalfEdgeIndices.size();
            }

            // Initialize attributes
            Attribute::init(mesh, attributeInitializerInfo);
        }

        Info extract(const SurfaceTriangularMesh& mesh) const {
            Info info;
            info.numVertices = mesh._vertices.size();
            const size_t numTriangles = mesh._triangles.size();
            info.triangleVertexIndexList.resize(numTriangles);

            for(size_t ti = 0; ti < numTriangles; ++ti) {
                size_t i = 0;
                mesh.forEachHalfEdgeInTriangle(ti, [&mesh, ti, &i, &info](size_t hei) {
                    info.triangleVertexIndexList[ti][i++] = mesh.target(hei);
                });
            }

            info.attributeInitializerInfo = Attribute::extract(mesh);

            return info;
        }
    };

    // Initialize the meshwork using triangle vertex index lists. Throws on error.
    template< typename Initializer, typename... Args > void init(Args&&... args) {
        _clear(); // Clear all the current topology
        Initializer().init(*this, std::forward<Args>(args)...);
    }
    template< typename Initializer > auto extract() const {
        return Initializer().extract(*this);
    }

    bool updateClosedness() {
        _isClosed = true;
        for(const auto& he : _halfEdges) if(!he.hasOpposite) { _isClosed = false; break; }
        return _isClosed;
    }
    bool isClosed()const noexcept { return _isClosed; }

    // Data accessors
    auto numVertices()  const noexcept { return _vertices.size(); }
    auto numHalfEdges() const noexcept { return _halfEdges.size(); }
    auto numEdges()     const noexcept { return _edges.size(); }
    auto numTriangles() const noexcept { return _triangles.size(); }
    const auto& getTriangles() const { return _triangles; }
    const auto& getHalfEdges() const { return _halfEdges; }
    const auto& getEdges()     const { return _edges; }
    const auto& getVertices()  const { return _vertices; }

    // Attribute accessor
    VertexAttribute&       getVertexAttribute(size_t index)       { return _vertices[index].attr; }
    const VertexAttribute& getVertexAttribute(size_t index) const { return _vertices[index].attr; }
    EdgeAttribute&       getEdgeAttribute(size_t index)       { return _edges[index].attr; }
    const EdgeAttribute& getEdgeAttribute(size_t index) const { return _edges[index].attr; }
    HalfEdgeAttribute&       getHalfEdgeAttribute(size_t index)       { return _halfEdges[index].attr; }
    const HalfEdgeAttribute& getHalfEdgeAttribute(size_t index) const { return _halfEdges[index].attr; }
    TriangleAttribute&       getTriangleAttribute(size_t index)       { return _triangles[index].attr; }
    const TriangleAttribute& getTriangleAttribute(size_t index) const { return _triangles[index].attr; }
    MetaAttribute&       getMetaAttribute()       { return _meta; }
    const MetaAttribute& getMetaAttribute() const { return _meta; }

    // Meshwork traverse
    bool hasOpposite(size_t halfEdgeIndex) const { return _halfEdges[halfEdgeIndex].hasOpposite; }
    size_t opposite(size_t halfEdgeIndex) const { return _halfEdges[halfEdgeIndex].oppositeHalfEdgeIndex; }
    size_t next(size_t halfEdgeIndex) const { return _halfEdges[halfEdgeIndex].nextHalfEdgeIndex; }
    size_t prev(size_t halfEdgeIndex) const { return _halfEdges[halfEdgeIndex].prevHalfEdgeIndex; }
    size_t triangle(size_t halfEdgeIndex) const { return _halfEdges[halfEdgeIndex].triangleIndex; }
    size_t target(size_t halfEdgeIndex) const { return _halfEdges[halfEdgeIndex].targetVertexIndex; }
    size_t edge(size_t halfEdgeIndex) const { return _halfEdges[halfEdgeIndex].edgeIndex; }

    size_t degree(size_t vertexIndex) const { return _vertices[vertexIndex].degree; }

    // Mesh neighbor iterators
    template< typename Func > void forEachHalfEdgeTargetingVertex(const Vertex& v, Func&& func) const {
        size_t hei0 = v.halfEdgeIndex;
        size_t hei = hei0;
        do {
            func(hei);
            if(!hasOpposite(hei)) break;
            hei = prev(opposite(hei));
        } while(hei != hei0);
    }
    template< typename Func > void forEachHalfEdgeTargetingVertex(size_t vi, Func&& func) const {
        forEachHalfEdgeTargetingVertex(_vertices[vi], std::forward<Func>(func));
    }
    template< typename Func > void forEachHalfEdgeInTriangle(const Triangle& t, Func&& func) const {
        size_t hei0 = t.halfEdgeIndex;
        size_t hei = hei0;
        do {
            func(hei);
            hei = next(hei);
        } while(hei != hei0);
    }
    template< typename Func > void forEachHalfEdgeInTriangle(size_t ti, Func&& func) const {
        forEachHalfEdgeInTriangle(_triangles[ti], std::forward<Func>(func));
    }
    template< typename Func > void forEachHalfEdgeInEdge(const Edge& e, Func&& func) const {
        size_t hei0 = e.halfEdgeIndex;
        func(hei0);
        if(hasOpposite(hei0)) func(opposite(hei0));
    }
    template< typename Func > void forEachHalfEdgeInEdge(size_t ei, Func&& func) const {
        forEachHalfEdgeInEdge(_edges[ei], std::forward<Func>(func));
    }

    // The following are basic mesh topology operators

    // Vertex insertion on edge.
    template< typename InsertionMethod >
    struct VertexInsertionOnEdge {
        static constexpr int deltaNumVertex = 1;

        template< typename AttributeSetter >
        void operator()(SurfaceTriangularMesh& mesh, size_t edgeIndex, AttributeSetter&& as)const {
            auto& edges = mesh._edges;
            auto& halfEdges = mesh._halfEdges;
            auto& vertices = mesh._vertices;
            auto& triangles = mesh._triangles;

            // Get index of current elements
            const size_t ohei       = edges[edgeIndex].halfEdgeIndex;
            const size_t ohei_n     = mesh.next(ohei);
            const size_t ohei_p     = mesh.prev(ohei);
            const size_t ohei_o     = mesh.opposite(ohei);
            const size_t ohei_on    = mesh.next(ohei_o);
            const size_t ohei_op    = mesh.prev(ohei_o);
            const size_t oti0       = mesh.triangle(ohei);
            const size_t oti2       = mesh.triangle(ohei_o);
            const size_t vi0        = mesh.target(ohei);
            const size_t vi1        = mesh.target(ohei_n);
            const size_t vi2        = mesh.target(ohei_o);
            const size_t vi3        = mesh.target(ohei_on);

            // Create new elements
            const InsertionMethod insertionMethod{ vi0, vi2 };
            const size_t vi     = mesh._newVertex(insertionMethod);
            const size_t ei2    = mesh._newEdge(insertionMethod); // New edge created by splitting
            const size_t hei0_o = mesh._newHalfEdge(insertionMethod); // Targeting new vertex, oppositing ohei
            const size_t hei2_o = mesh._newHalfEdge(insertionMethod); // Targeting new vertex, oppositing ohei_o
            const size_t ei1    = mesh._newEdge(insertionMethod); // New edge cutting t0
            const size_t hei1   = mesh._newHalfEdge(insertionMethod); // Leaving new vertex
            const size_t hei1_o = mesh._newHalfEdge(insertionMethod); // Targeting new vertex
            const size_t ei3    = mesh._newEdge(insertionMethod); // New edge cutting t2
            const size_t hei3   = mesh._newHalfEdge(insertionMethod); // Leaving new vertex
            const size_t hei3_o = mesh._newHalfEdge(insertionMethod); // Targeting new vertex
            const size_t ti1    = mesh._newTriangle(insertionMethod);
            const size_t ti3    = mesh._newTriangle(insertionMethod);

            // Adjust vertex
            vertices[vi].halfEdgeIndex = hei0_o;
            halfEdges[hei0_o].targetVertexIndex = vi;
            halfEdges[hei2_o].targetVertexIndex = vi;
            halfEdges[hei1_o].targetVertexIndex = vi;
            halfEdges[hei3_o].targetVertexIndex = vi;
            halfEdges[hei1].targetVertexIndex = vi1;
            halfEdges[hei3].targetVertexIndex = vi3;

            vertices[vi].degree = 4;
            ++vertices[vi1].degree;
            ++vertices[vi3].degree;

            // Adjust triangle
            mesh._registerTriangle(oti0, ohei,   ohei_n,  hei1_o);
            mesh._registerTriangle(ti1,  hei1,   ohei_p,  hei2_o);
            mesh._registerTriangle(oti2, ohei_o, ohei_on, hei3_o);
            mesh._registerTriangle(ti3,  hei3,   ohei_op, hei0_o);

            // Adjust edge
            mesh._registerEdge(edgeIndex, ohei,   hei0_o);
            mesh._registerEdge(ei1,       hei1,   hei1_o);
            mesh._registerEdge(ei2,       ohei_o, hei2_o);
            mesh._registerEdge(ei3,       hei3,   hei3_o);

            // Update attributes of affected elements
            as(
                mesh,
                {oti0, ti1, oti2, ti3},
                {vi, vi0, vi1, vi2, vi3},
                {edgeIndex, ei1, ei2, ei3}
            );
        }

        void operator()(SurfaceTriangularMesh& mesh, size_t edgeIndex)const {
            this->operator()(mesh, edgeIndex, [](
                SurfaceTriangularMesh& mesh,
                std::array<size_t, 4> tis,
                std::array<size_t, 5> vis,
                std::array<size_t, 4> eis
            ) {});
        }
    };
    // Edge collapse
    struct EdgeCollapse {
        static constexpr int deltaNumVertex = -1;

        // The target of the halfedge ohei will be preserved
        // Notice that halfedge index (not edge index) is used in this function.
        template< typename AttributeSetter >
        void operator()(SurfaceTriangularMesh& mesh, size_t ohei, AttributeSetter&& as)const {
            auto& edges = mesh._edges;
            auto& halfEdges = mesh._halfEdges;
            auto& vertices = mesh._vertices;
            auto& triangles = mesh._triangles;

            // Preconditions should be handled by the caller

            // Get index of current elements
            const size_t oei = mesh.edge(ohei);
            const size_t ohei_n = mesh.next(ohei);
            const size_t ohei_p = mesh.prev(ohei);
            const size_t ohei_o = mesh.opposite(ohei);
            const size_t ohei_on = mesh.next(ohei_o);
            const size_t ohei_op = mesh.prev(ohei_o);
            const size_t ot0 = mesh.triangle(ohei);
            const size_t ot1 = mesh.triangle(ohei_o);
            const size_t ov0 = mesh.target(ohei); // Will collapse to this vertex
            const size_t ov1 = mesh.target(ohei_o); // Will be removed
            const size_t oei1 = mesh.edge(ohei_n); // Will collapse to this edge
            const size_t oei2 = mesh.edge(ohei_p); // Will be removed
            const size_t oei3 = mesh.edge(ohei_op); // Will collapse on this edge
            const size_t oei4 = mesh.edge(ohei_on); // Will be removed

            // Retarget all halfedges pointing v1 to v0
            const size_t hei_begin = mesh.opposite(ohei_on); // Changed halfedge begin
            const size_t hei_end = mesh.opposite(ohei_n);    // Changed halfedge end
            for(size_t hei1 = hei_begin; hei1 != ohei_p; hei1 = mesh.opposite(mesh.next(hei1))) {
                halfEdges[hei1].targetVertexIndex = ov0;
            }
            vertices[ov0].halfEdgeIndex = hei_begin;

            // Collapse edges
            mesh._registerEdge(oei1, mesh.opposite(ohei_n), mesh.opposite(ohei_p));
            mesh._registerEdge(oei3, mesh.opposite(ohei_op), mesh.opposite(ohei_on));

            // Adjust vertex degrees
            vertices[ov0].degree += vertices[ov1].degree - 4;
            --vertices[mesh.target(ohei_n)].degree;
            --vertices[mesh.target(ohei_on)].degree;

            // Remove elements
            mesh._removeVertex(ov1);
            mesh._removeEdge(oei);
            mesh._removeEdge(oei2);
            mesh._removeEdge(oei4);
            mesh._removeHalfEdge(ohei);   mesh._removeHalfEdge(ohei_n);  mesh._removeHalfEdge(ohei_p);
            mesh._removeHalfEdge(ohei_o); mesh._removeHalfEdge(ohei_on); mesh._removeHalfEdge(ohei_op);
            mesh._removeTriangle(ot0);
            mesh._removeTriangle(ot1);

            // Update attributes for affected elements
            as(
                mesh,
                hei_begin, hei_end, // Range of changed halfedges targeting ov0, ordered clockwise
                ov0
            );
        }

        void operator()(SurfaceTriangularMesh& mesh, size_t ohei)const {
            this->operator()(mesh, ohei, [](
                SurfaceTriangularMesh& mesh,
                size_t hei_begin, size_t hei_end,
                size_t ov0
            ) {});
        }
    };
    // Edge flip
    struct EdgeFlip {
        static constexpr int deltaNumVertex = 0;

        template< typename AttributeSetter >
        void operator()(SurfaceTriangularMesh& mesh, size_t edgeIndex, AttributeSetter&& as) const {
            auto& edges = mesh._edges;
            auto& halfEdges = mesh._halfEdges;
            auto& vertices = mesh._vertices;
            auto& triangles = mesh._triangles;

            // Preconditions (like topology) should be handled by the caller.

            // Get index of current elements
            const size_t ohei = edges[edgeIndex].halfEdgeIndex;
            const size_t ohei_n = mesh.next(ohei);
            const size_t ohei_p = mesh.prev(ohei);
            const size_t ohei_o = mesh.opposite(ohei);
            const size_t ohei_on = mesh.next(ohei_o);
            const size_t ohei_op = mesh.prev(ohei_o);
            const size_t ov0 = mesh.target(ohei);
            const size_t ov1 = mesh.target(ohei_n);
            const size_t ov2 = mesh.target(ohei_o);
            const size_t ov3 = mesh.target(ohei_on);
            const size_t ot0 = mesh.triangle(ohei);
            const size_t ot1 = mesh.triangle(ohei_o);

            // Retarget vertices
            halfEdges[ohei].targetVertexIndex = ov1;
            halfEdges[ohei_o].targetVertexIndex = ov3;
            vertices[ov0].halfEdgeIndex = ohei_op;
            vertices[ov2].halfEdgeIndex = ohei_p;

            --vertices[ov0].degree;
            --vertices[ov2].degree;
            ++vertices[ov1].degree;
            ++vertices[ov3].degree;

            // Remake triangles
            mesh._registerTriangle(ot0, ohei, ohei_p, ohei_on);
            mesh._registerTriangle(ot1, ohei_o, ohei_op, ohei_n);

            // Update attributes of affected elements
            as(
                mesh,
                {ot0, ot1},
                {ov0, ov1, ov2, ov3}
            );

        }

        void operator()(SurfaceTriangularMesh& mesh, size_t edgeIndex) const {
            this->operator()(mesh, edgeIndex, [](
                SurfaceTriangularMesh& mesh,
                std::array<size_t, 2> tis,
                std::array<size_t, 4> vis
            ) {});
        }

    };

    // triangle subdivision, introduces 3 new vertices
    template< typename InsertionMethod > struct [[deprecated]] TriangleSubdivision {
        static constexpr int deltaNumVertex = 3;
        void operator()(SurfaceTriangularMesh& mesh, size_t triangleIndex) const {
            auto& edges = mesh._edges;
            auto& halfEdges = mesh._halfEdges;
            auto& vertices = mesh._vertices;
            auto& triangles = mesh._triangles;

            // TODO pre-conditions

            const size_t ohei = triangles[triangleIndex].halfEdgeIndex;
            const size_t ohei_n = mesh.next(ohei);
            const size_t ohei_p = mesh.prev(ohei);
            VertexInsertionOnEdge<InsertionMethod>()(mesh, mesh.edge(ohei)); // ohei, ohei_n, ohei_p still have the same target

            // Get the edge for future flipping.
            const size_t eFlip = mesh.edge(mesh.prev(ohei));

            VertexInsertionOnEdge<InsertionMethod>()(mesh, mesh.edge(ohei_n));
            VertexInsertionOnEdge<InsertionMethod>()(mesh, mesh.edge(ohei_p));
            EdgeFlip()(mesh, eFlip);
        }
    };

};

#endif
