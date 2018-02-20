#ifndef MEDYAN_MembraneRegion_h
#define MEDYAN_MembraneRegion_h

#include <array>
#include <vector>

// Forward Declarations
class MembraneHierarchy;
class Boundary;

/******************************************************************************
 * 
 * Membrane Region is a class for management of regions enclosed by membranes
 * and boundaries.
 * 
******************************************************************************/

class MembraneRegion {
private:
    Boundary* _boundary = nullptr; ///< The boundary of the playground
    std::vector<MembraneHierarchy*> _hierOut; ///< Their union is the outer limit
    std::vector<MembraneHierarchy*> _hierIn;  ///< Their union is the inner limit

public:
    /// Construct the region with one membrane
    MembraneRegion(MembraneHierarchy* hier, bool excludeChildren=false);

    /// Construct the region with the boundary
    MembraneRegion(Boundary* b): _boundary(b) {}
    MembraneRegion(Boundary* b, MembraneHierarchy* parentOfExcluded);

    /// Is point inside region
    bool contains(const std::array<double, 3>& point)const;
};

#endif