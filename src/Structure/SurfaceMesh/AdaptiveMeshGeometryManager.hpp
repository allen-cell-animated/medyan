#ifndef MEDYAN_Structure_SurfaceMesh_AdaptiveMeshGeometryManager_hpp
#define MEDYAN_Structure_SurfaceMesh_AdaptiveMeshGeometryManager_hpp

namespace adaptive_mesh {

template< typename Mesh > struct GeometryManager {
    static void computeAllTriangleNormals(Mesh& mesh) {
        const size_t numTriangles = mesh.numTriangles();

        for(typename Mesh::TriangleIndex ti { 0 }; ti < numTriangles; ++ti) {
            medyan::adaptiveComputeTriangleNormal(mesh, ti);
        }
    }

    static void computeAllAngles(Mesh& mesh) {
        const size_t numHalfEdges = mesh.numHalfEdges();

        for(typename Mesh::HalfEdgeIndex hei { 0 }; hei < numHalfEdges; ++hei) {
            medyan::adaptiveComputeAngle(mesh, hei);
        }
    }

    // Requires
    //   - Unit normals in triangles (geometric)
    //   - Angles in halfedges (geometric)
    static void computeAllVertexNormals(Mesh& mesh) {
        const size_t numVertices = mesh.numVertices();

        for(typename Mesh::VertexIndex vi {0}; vi < numVertices; ++vi) {
            medyan::adaptiveComputeVertexNormal(mesh, vi);
        }
    }
}; // End GeometryManager

} // namespace adaptive_mesh

#endif
