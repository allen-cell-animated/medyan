#ifndef MEDYAN_Structure_SurfaceMesh_MembraneMeshVertexSystem_hpp
#define MEDYAN_Structure_SurfaceMesh_MembraneMeshVertexSystem_hpp

namespace medyan {

// If vertices represent fixed coordinates in a material coordinate system,
// a vertex represents the location of a specific lipid molecule.
// The local area elasticity is usually a necessity in this case.
// An exception is the vertices on the border connected to a lipid
// reservoir, where the vertices stand for the boundary locations which are
// usually fixed in the ambient space.
//
// If vertices represent fixed coordinates in a normal coordinate system,
// a vertex is only a representative point on the surface, where the motion
// of the vertex must be in the local normal direction.
// Local area elasticity cannot be directly defined on mesh elements.
//
// If vertices represent fixed coordinates in a normal coordinate system,
// then the system is similar to the "normal" coordinate case, but without
// the requirement for vertices to move only in the normal directions.
enum class MembraneMeshVertexSystem { material, normal, general };

} // namespace medyan

#endif
