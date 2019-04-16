#ifndef MEDYAN_Structure_SurfaceMesh_AdaptiveMeshGeometryManager_hpp
#define MEDYAN_Structure_SurfaceMesh_AdaptiveMeshGeometryManager_hpp

namespace adaptive_mesh {

template< typename Mesh > struct GeometryManager {
    static void computeAllTriangleNormals(Mesh& mesh) {
        const size_t numTriangles = mesh.getTriangles().size();

        for(size_t ti = 0; ti < numTriangles; ++ti) {
            Mesh::AttributeType::adaptiveComputeTriangleNormal(mesh, ti);
        }
    }

    static void computeAllAngles(Mesh& mesh) {
        const size_t numHalfEdges = mesh.getHalfEdges().size();

        for(size_t hei = 0; hei < numHalfEdges; ++hei) {
            Mesh::AttributeType::adaptiveComputeAngle(mesh, hei);
        }
    }

    // Requires
    //   - Unit normals in triangles (geometric)
    //   - Angles in halfedges (geometric)
    static void computeAllVertexNormals(Mesh& mesh) {
        const size_t numVertices = mesh.getVertices().size();

        for(size_t vi = 0; vi < numVertices; ++vi) if(!mesh.isVertexOnBorder(vi)) {
            Mesh::AttributeType::adaptiveComputeVertexNormal(mesh, vi);
        }
    }
}; // End GeometryManager

} // namespace adaptive_mesh

#endif
