/*

Adaptive mesh algorithm

Implementation based on
"An Adaptive Mesh Algorithm for Evolving Surfaces: Simulations of Drop Breakup and Coalescence"
by Vittorio Cristini, Jerzy Blawzdziewicz and Michael Loewenberg.

Performs mesh relaxation and topological transformation.

Ending criteria:
    1. Relaxation terminates.
    2. in a local region, number of nodes n is near n_0.
    3. no more node reconnection.

The following values are stored and updated
    - L computed on each vertex.
        - L_0 is generally universal, as an upper limit
        - L_1 uses local curvature information.
        - L_2 and above are not used currently.
    - rho computed on each vertex. rho = c_0 / (alpha * L)^2
        - c_0 is a geometric factor (2 / sqrt(3))
        - alpha_0 is resolution of the surface (ref: 0.2-0.3)
    - rho_avg computed on each vertex. (1-ring weighted avg of rho)
    - A_0 computed an each triangle. Uses rho_avg on vertices.

For mesh relaxation, the following values are also needed
    - l_0' computed on each vertex. Uses weighted sum of neighboring A_0.
    - l_0 computed on each edge. An average of l_0' on the two vertex.
    - normal computed on vertices.
    - v (velocity) on vertices, Uses l_0 and n.

For topological transformation, the following values are also needed
    - S_loc computed on each vertex
    - S_0 computed on a local region. Uses A_0 in the region.
    - n_0 computed on a local region. n_0 = n * S_loc / S_0. (fraction of border nodes will be counted to n).

We might need to redistribute the surface area globally or locally depending on the operation,
if we are using patched area energy computation (i.e. use sum of local area energy instead of
energy from sum of area)

The algorithm was not introduced explicity in the article, and we'll formalize it as follows

Init: Find per-element quantities for all elements.
Loop
    Global relaxation
    Update per-element quantities and local averaged quantities
    For all places that does not satisfy criteria 2, 3
        Local topological transformation
        Local relaxation
        Update affected per-element quantities and local averaged quantites
    End
Until all criteria are met

It would be easier if the implementation and the mesh representation are coupled with
the current surface meshwork system. But additional variables should be introduced and
it might not be appropriate to mess them up with the original structure.

Overlaying a new set of variables with the current implementation of the meshwork might
be a good idea as well.

*/

#ifndef ADAPTIVE_MESH_HPP
#define ADAPTIVE_MESH_HPP

#include <vector>

using Float = double;

struct MeshForce {
    const std::vector< Vertex* >& vertices;
    std::vector< char > mask;
    std::vector< size_t > activeVertexIndices;
    some_array operator()(some_array coord) {
        // The size of coord must be exactly 3 times the size of vertex array.
        for(size_t i : activeVertexIndices) {
            Vector3 p = vertices[i] -> point;
            for(Vertex* v : vertices[i] -> neighborVertices) {
                if
            }
        }
    }
};

// 2nd order Runge Kutta method on a set of coordinates.
void rk2(coord, calc_force, Float dt, Float epsilon) {
    force = calc_force(coord);
    coord_halfway = coord + force * dt / 2;
    force_halfway = calc_force(coord_halfway);
    coord += force_halfway * dt;
}

#endif
