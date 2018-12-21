#ifndef MEDYAN_SurfaceMesh_hpp
#define MEDYAN_SurfaceMesh_hpp

#include <algorithm>
#include <array>
#include <cstddef> // ptrdiff_t
#include <iterator>
#include <set>
#include <type_traits>
#include <vector>

/******************************************************************************
The data structure for an orientable, manifold 2d triangular meshwork in 3d
space.

The connectivity is similar to CGAL library, which is halfedge based.

target:     The vertex that the halfedge points to.
opposite:   The halfedge on the same edge but in the opposite direction.
face:       The triangle that the halfedge is circling ccw.
next:       The next halfedge in the current face circling ccw.
prev:       The previous halfedge in the current face circling ccw.
edge:       The undirected edge of this halfedge.

All the other elements must have at least one index pointing to a halfedge.
******************************************************************************/

template< typename T > class DeletableVector {
public:
    struct Deletable {
        bool markedAsDeleted = false;
        T data;
    };
    struct IndexMove {
        size_t from, to;
        bool valid;
    };

private:
    std::vector< Deletable > _value;
    std::set< size_t > _deletedIndices;

    template< bool constIt > class IndexIteratorImpl {
    public:
        using iterator_category = std::bidirectional_iterator_tag;
        using value_type = int;
        using difference_type = std::ptrdiff_t;
        using pointer = typename std::conditional_t< constIt, const int*, int* >;
        using reference = typename std::conditional_t< constIt, const int&, int& >;

    private:
        int _pos; // index in _value
        std::vector< Deletable >* _storage;

    public:
        explicit IndexIteratorImpl(int pos, std::vector< Deletable >* storage) : _pos(pos), _storage(storage) {
            while(_pos >= 0 && _pos < _storage->size() && (*_storage)[_pos].markedAsDeleted) ++_pos;
        }

        reference operator*()const { return _pos; }
        pointer operator->()const { return &_pos; }

        IndexIteratorImpl& operator++() {
            do ++_pos; while(_pos >= 0 && _pos < _storage->size() && (*_storage)[_pos].markedAsDeleted);
            return *this;
        }
        IndexIteratorImpl operator++(int) {
            IndexIteratorImpl res = *this;
            do ++_pos; while(_pos >= 0 && _pos < _storage->size() && (*_storage)[_pos].markedAsDeleted);
            return res;
        }
        IndexIteratorImpl& operator--() {
            do --_pos; while(_pos >= 0 && _pos < _storage->size() && (*_storage)[_pos].markedAsDeleted);
            return *this;
        }
        IndexIteratorImpl operator--(int) {
            IndexIteratorImpl res = *this;
            do --_pos; while(_pos >= 0 && _pos < _storage->size() && (*_storage)[_pos].markedAsDeleted);
            return res;
        }

        template< bool rhsConstIt >
        bool operator==(const IndexIteratorImpl< rhsConstIt >& rhs)const { return _storage == rhs._storage && _pos == rhs.pos; }
    };

public:

    using iterator       = IndexIteratorImpl< false >;
    using const_iterator = IndexIteratorImpl< true >;
    iterator       begin() noexcept       { return iterator(0, _value); }
    const_iterator begin() const noexcept { return const_iterator(0, _value); }
    iterator       end() noexcept       { return iterator(_value.size(), _value); }
    const_iterator end() const noexcept { return const_iterator(_value.size(), _value); }

    // Insert a new element. Returns the new index.
    size_t insert() {
        if(_deletedIndices.empty()) {
            _value.emplace_back();
            return _value.size() - 1;
        } else {
            auto it = _deletedIndices.cbegin(); // minimum
            const size_t res = *it;
            _value[res].markedAsDeleted = false;
            _deletedIndices.erase(it);
            return res;
        }
    }

    void erase(size_t index) {
        _value[index].markedAsDeleted = true;
        _deletedIndices.insert(index);
    }

    // Make the indices an exact copy from another DeletableVector
    template< typename U >
    void resizeFrom(const DeletableVector<U>& rhs) {
        _deletedIndices = rhs._deletedIndices;
        _value.resize(rhs._value.size());
        for(size_t i : _deletedIndices) _value[i]._markedAsDeleted = true;
    }

    bool isDeleted(size_t index) const { return _value[index].markedAsDeleted; }

    // Move the last undeleted element to the first deleted location.
    // Returns the index of the element moved.
    IndexMove makeMove() {
        IndexMove res;

        // Remove trailing deleted elements.
        while(!_value.empty() && _value.back().markedAsDeleted) {
            _value.pop_back();
            _deletedIndices.erase(_deletedIndices.crbegin()); // Remove maximum index
        }

        if(_deletedIndices.empty()) {
            res.valid = false; // Invalid move
        } else {
            const auto it = _deletedIndices.cbegin(); // minimum index
            res.from = _value.size() - 1;
            res.to = *it;
            res.valid = true;
            _value[res.to] = std::move(_value[res.from]);
            _value.pop_back();
            _deletedIndices.erase(it);
        }

        return res;
    }

    size_t size_raw() const noexcept { return _value.size(); }
    size_t size() const { return _value.size() - _deletedIndices.size(); }

    T&       operator[](size_t index)       { return _value[index].data; }
    const T& operator[](size_t index) const { return _value[index].data; }

    // This function should only be called during initialization
    std::vector< Deletable >& getValue() { return _value; }

};

template< typename Attribute > class SurfaceTriangularMeshBase {
public:

    using VertexAttribute   = typename Attribute::VertexAttribute;
    using EdgeAttribute     = typename Attribute::EdgeAttribute;
    using HalfEdgeAttribute = typename Attribute::HalfEdgeAttribute;
    using TriangleAttribute = typename Attribute::TriangleAttribute;
    using MetaAttribute     = typename Attribute::MetaAttribute;

    // The elements should be trivially copyable.
    struct Vertex {
        size_t halfEdgeIndex; // Only one HalfEdge targeting the vertex is needed.
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

protected:

    DeletableVector<Triangle> _triangles; // collection of triangles
    DeletableVector<HalfEdge> _halfEdges; // collection of halfedges
    DeletableVector<Edge>     _edges;     // collection of edges
    DeletableVector<Vertex>   _vertices;  // collection of vertices

    MetaAttribute _meta;

    bool _isClosed = true; // Whether the meshwork is topologically closed
    int _genus = 0; // Genus of the surface. Normally 0, as for a topologically spherical shape

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
        Attribute::newVertex(_meta, *this, index, op);
        return index;
    }

public:

    // Constructors
    SurfaceTriangularMeshBase(const MetaAttribute& meta) : _meta(meta) {}

    // Initialize the meshwork using triangle vertex index lists. Throws on error.
    void init(
        size_t numVertices,
        const std::vector< std::array< size_t, 3 > >& triangleVertexIndexList
    ) {
        _vertices.getValue().resize(numVertices);
        const size_t numTriangles = triangleVertexIndexList.size();
        _triangles.getValue().resize(numTriangles);
        const size_t numHalfEdges = 3 * numTriangles;
        _halfEdges.getValue().reserve(numHalfEdges);
        _edges.getValue().reserve(numHalfEdges / 2); // Might be more than this number with borders.

        struct VertexAdditionalInfo {
            bool hasTargetingHalfEdge = false;
            std::vector< size_t > leavingHalfEdgeIndices;
        };
        std::vector< VertexAdditionalInfo > vai(numVertices);

        for(size_t ti = 0; ti < numTriangles; ++ti) {
            const auto& t = triangleVertexIndexList[ti];
            _triangles[ti].halfEdgeIndex = _halfEdges.size_raw(); // The next inserted halfedge index
            for(size_t i = 0; i < 3; ++i) {
                const size_t hei = _halfEdges.insert();
                HalfEdge& he = _halfEdges[hei];
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
                    [this, leftVertexIndex](size_t leavingHalfEdgeIndex) {
                        return leftVertexIndex == _halfEdges[leavingHalfEdgeIndex].targetVertexIndex;
                    }
                );
                if(findRes == vai[t[i]].leavingHalfEdgeIndices.end()) {
                    // opposite not found
                    _edges[_edges.insert()].halfEdgeIndex = hei;
                    he.edgeIndex = _edges.size_raw() - 1;
                } else {
                    // opposite found
                    he.hasOpposite = true;
                    he.oppositeHalfEdgeIndex = *findRes;
                    he.edgeIndex = _halfEdges[he.oppositeHalfEdgeIndex].edgeIndex;

                    _halfEdges[he.oppositeHalfEdgeIndex].hasOpposite = true;
                    _halfEdges[he.oppositeHalfEdgeIndex].oppositeHalfEdgeIndex = hei;
                }

                // Set vertex half edge index
                if(!vai[t[i]].hasTargetingHalfEdge) {
                    _vertices[t[i]].halfEdgeIndex = hei;
                    vai[t[i]].hasTargetingHalfEdge = true;
                }
            }
        }

        // TODO: validation
    }

    // Attribute accessor
    VertexAttribute&       getVertexAttribute(size_t index)       { return _vertices[index].attr; }
    const VertexAttribute& getVertexAttribute(size_t index) const { return _vertices[index].attr; }
    EdgeAttribute&       getEdgeAttribute(size_t index)       { return _edges[index].attr; }
    const EdgeAttribute& getEdgeAttribute(size_t index) const { return _edges[index].attr; }
    HalfEdgeAttribute&       getHalfEdgeAttribute(size_t index)       { return _halfEdges[index].attr; }
    const HalfEdgeAttribute& getHalfEdgeAttribute(size_t index) const { return _halfEdges[index].attr; }
    TriangleAttribute&       getTriangleAttribute(size_t index)       { return _triangles[index].attr; }
    const TriangleAttribute& getTriangleAttribute(size_t index) const { return _triangles[index].attr; }

    // Meshwork traverse
    bool hasOpposite(size_t halfEdgeIndex) const { return _halfEdges[halfEdgeIndex].hasOpposite; }
    size_t opposite(size_t halfEdgeIndex) const { return _halfEdges[halfEdgeIndex].oppositeHalfEdgeIndex; }
    size_t next(size_t halfEdgeIndex) const { return _halfEdges[halfEdgeIndex].nextHalfEdgeIndex; }
    size_t prev(size_t halfEdgeIndex) const { return _halfEdges[halfEdgeIndex].prevHalfEdgeIndex; }
    size_t triangle(size_t halfEdgeIndex) const { return _halfEdges[halfEdgeIndex].triangleIndex; }
    size_t target(size_t halfEdgeIndex) const { return _halfEdges[halfEdgeIndex].targetVertexIndex; }
    size_t edge(size_t halfEdgeIndex) const { return _halfEdges[halfEdgeIndex].edgeIndex; }

    // The following are basic mesh topology operators

    // Vertex insertion on edge.
    struct VertexInsertionOnEdge {
        static constexpr int deltaNumVertex = 1;

        struct InsertMid { size_t v0, v1; };

        void operator()(SurfaceTriangularMeshBase& mesh, size_t edgeIndex)const {
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
            const size_t vi     = mesh._newVertex(InsertMid{vi0, vi2});
            const size_t ei2    = edges.insert(); // New edge created by splitting
            const size_t hei0_o = halfEdges.insert(); // Targeting new vertex, oppositing ohei
            const size_t hei2_o = halfEdges.insert(); // Targeting new vertex, oppositing ohei_o
            const size_t ei1    = edges.insert(); // New edge cutting t0
            const size_t hei1   = halfEdges.insert(); // Leaving new vertex
            const size_t hei1_o = halfEdges.insert(); // Targeting new vertex
            const size_t ei3    = edges.insert(); // New edge cutting t2
            const size_t hei3   = halfEdges.insert(); // Leaving new vertex
            const size_t hei3_o = halfEdges.insert(); // Targeting new vertex
            const size_t ti1    = triangles.insert();
            const size_t ti3    = triangles.insert();

            // Adjust vertex
            vertices[vi].halfEdgeIndex = hei0_o;
            halfEdges[hei0_o].targetVertexIndex = vi;
            halfEdges[hei2_o].targetVertexIndex = vi;
            halfEdges[hei1_o].targetVertexIndex = vi;
            halfEdges[hei3_o].targetVertexIndex = vi;
            halfEdges[hei1].targetVertexIndex = vi1;
            halfEdges[hei3].targetVertexIndex = vi3;

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
        }
    };
    // Edge collapse
    struct EdgeCollapse {
        static constexpr int deltaNumVertex = -1;
        void operator()(SurfaceTriangularMeshBase& mesh, size_t edgeIndex)const {
            auto& edges = mesh._edges;
            auto& halfEdges = mesh._halfEdges;
            auto& vertices = mesh._vertices;
            auto& triangles = mesh._triangles;

            // TODO preconditions

            // Get index of current elements
            const size_t ohei = edges[edgeIndex].halfEdgeIndex;
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
            for(size_t hei1 = mesh.opposite(ohei_on); hei1 != ohei_p; hei1 = mesh.opposite(mesh.next(hei1))) {
                halfEdges[hei1].targetVertexIndex = ov0;
            }
            vertices[ov0].halfEdgeIndex = ohei_on;

            // Collapse edges
            mesh._registerEdge(oei1, mesh.opposite(ohei_n), mesh.opposite(ohei_p));
            mesh._registerEdge(oei3, mesh.opposite(ohei_op), mesh.opposite(ohei_on));

            // Remove elements
            vertices.erase(ov1);
            edges.erase(edgeIndex);
            edges.erase(oei2);
            edges.erase(oei4);
            halfEdges.erase(ohei);   halfEdges.erase(ohei_n);  halfEdges.erase(ohei_p);
            halfEdges.erase(ohei_o); halfEdges.erase(ohei_on); halfEdges.erase(ohei_op);
            triangles.erase(ot0);
            triangles.erase(ot1);
        }
    };
    // Edge flip
    struct EdgeFlip {
        static constexpr int deltaNumVertex = 0;
        void operator()(SurfaceTriangularMeshBase& mesh, size_t edgeIndex) {
            auto& edges = mesh._edges;
            auto& halfEdges = mesh._halfEdges;
            auto& vertices = mesh._vertices;
            auto& triangles = mesh._triangles;

            // TODO precondition

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

            // Remake triangles
            mesh._registerTriangle(ot0, ohei, ohei_p, ohei_on);
            mesh._registerTriangle(ot1, ohei_o, ohei_op, ohei_n);

        }
    };

    // triangle subdivision, introduces 3 new vertices
    struct TriangleSubdivision {
        static constexpr int deltaNumVertex = 3;
        void operator()(SurfaceTriangularMeshBase& mesh, size_t triangleIndex) {
            auto& edges = mesh._edges;
            auto& halfEdges = mesh._halfEdges;
            auto& vertices = mesh._vertices;
            auto& triangles = mesh._triangles;

            // TODO pre-conditions

            const size_t ohei = triangles[triangleIndex].halfEdgeIndex;
            const size_t ohei_n = mesh.next(ohei);
            const size_t ohei_p = mesh.prev(ohei);
            VertexInsertionOnEdge()(mesh, mesh.edge(ohei)); // ohei, ohei_n, ohei_p still have the same target

            // Get the edge for future flipping.
            const size_t eFlip = mesh.edge(mesh.prev(ohei));

            VertexInsertionOnEdge()(mesh, mesh.edge(ohei_n));
            VertexInsertionOnEdge()(mesh, mesh.edge(ohei_p));
            EdgeFlip()(mesh, eFlip);
        }
    };

};

// Type GeometricAttribute::VertexData must implement getCoordinate() function
template< typename GeometricAttribute > class SurfaceTriangularMesh : public SurfaceTriangularMeshBase< GeometricAttribute > {
public:
    using Base = SurfaceTriangularMeshBase< GeometricAttribute >;

    using VertexAttribute   = typename Base::VertexAttribute;
    using EdgeAttribute     = typename Base::EdgeAttribute;
    using HalfEdgeAttribute = typename Base::HalfEdgeAttribute;
    using TriangleAttribute = typename Base::TriangleAttribute;
    using MetaAttribute     = typename Base::MetaAttribute;

    struct GeometricVertexInit {};

    using coordinate_type = typename GeometricAttribute::VertexAttribute::coordinate_type;

    // Constructors
    SurfaceTriangularMesh(const MetaAttribute& meta) : Base(meta) {}

    void init(
        const std::vector< coordinate_type >& vertexCoordinateList,
        const std::vector< std::array< size_t, 3 > >& triangleVertexIndexList
    ) {
        const size_t numVertices = vertexCoordinateList.size();

        Base::init(numVertices, triangleVertexIndexList);

        for(size_t i = 0; i < numVertices; ++i) Attribute::newVertex(_meta, *this, i, vertexCoordinateList[i], GeometricVertexInit{});
        for(auto& ei : _edges) Attribute::newEdge(_meta, *this, ei, GeometricVertexInit{});
        // TODO Perform other initialization steps


    }

};



#endif
