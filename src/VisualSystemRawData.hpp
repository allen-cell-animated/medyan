#ifndef MEDYAN_VisualSystemRawData_hpp
#define MEDYAN_VisualSystemRawData_hpp

#include <array>
#include <cstdint>
#include <mutex>
#include <vector>

#include "Structure/Bead.h"
#include "Util/Math/Vec.hpp"

namespace visual {

namespace raw_data_cat {

    using Type = std::uint_fast8_t;

    constexpr Type none            = 0;
    constexpr Type beadPosition    = 1 << 0;
    constexpr Type beadConnection  = 1 << 1;
    constexpr Type compartment     = 1 << 2;
    constexpr Type concentration   = 1 << 3;

} // namespace raw_data_cat

// SystemRawData is a stateful data container, used for storing system data
// that is useful for visualization.
//
// The raw data is updated in the simulation thread with system data.
// The data can be consumed so that the update state is reset.
struct SystemRawData {
    using V3 = mathfunc::Vec3;

    struct MembraneData {
        std::vector< V3 > vertexCoords;
        std::vector< std::array< size_t, 3 > > triangleVertexIndices; // Index for this membrane only (start from 0)
    };
    struct FilamentData {
        std::vector< V3 > beadCoords;
    };

    // Synchronization and states
    std::mutex me;

    raw_data_cat::Type updated = raw_data_cat::none; // Should be reset by the consumer

    // Data
    mathfunc::Vec< 3, size_t > compartmentNum;
    mathfunc::Vec3             compartmentSize;

    std::vector< MembraneData > membraneData;
    std::vector< FilamentData > filamentData; // [Filament Idx][Bead position in filament]

    std::vector< std::array< mathfunc::Vec3, 2 > > linkerCoords;
    std::vector< std::array< mathfunc::Vec3, 2 > > motorCoords;
    std::vector< std::array< mathfunc::Vec3, 2 > > brancherCoords;

};

// Function to copy the system data to raw data
//
// Note:
//   - This function should only be called on the simulation thread
//
// Input:
//   - data: the raw data object
//   - updated: the part of the system that's needs to be updated
//   - ignoreDataInUse: if set to true, copy would be skipped if the data is
//     being used by other threads
//
// Returns whether the copy is actually successfully made.
bool copySystemData(
    SystemRawData& data,
    raw_data_cat::Type updated,
    bool ignoreDataInUse = true
);

} // namespace visual

#endif
