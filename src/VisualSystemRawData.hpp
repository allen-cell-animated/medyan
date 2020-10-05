#ifndef MEDYAN_VisualSystemRawData_hpp
#define MEDYAN_VisualSystemRawData_hpp

#include <array>
#include <cstdint>
#include <mutex>
#include <vector>

#include "Structure/Bead.h"
#include "Util/Math/Vec.hpp"
#include "Visual/FrameData.hpp"
#include "Visual/VisualElement.hpp"

namespace medyan::visual {

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

    // Synchronization and states
    std::mutex me;

    raw_data_cat::Type updated = raw_data_cat::none; // Should be reset by the consumer

    // meta data
    DisplayTypeMap             displayTypeMap;

    // Data
    DisplayFrame               frameData;

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


// Note:
//   - This function should only be called on the visual thread.
inline void updateRealtimeMeshData(
    std::vector< ProfileWithMeshData >& profileData,
    SystemRawData&                      rawData
) {
    std::lock_guard< std::mutex > rawGuard(rawData.me);

    const auto updateMeshData = [&](MeshData& meshData) {
        return Overloaded {
            [&](const MembraneProfile& profile) {
                if(rawData.updated & raw_data_cat::beadPosition) {
                    meshData = createMembraneMeshData(
                        select(rawData.frameData.membranes, profile.selector),
                        profile.displaySettings
                    );
                }
            },
            [&](const FilamentProfile& profile) {
                if(rawData.updated & raw_data_cat::beadPosition) {
                    meshData = createFilamentMeshData(
                        select(rawData.frameData.filaments, profile.selector),
                        profile.displaySettings
                    );
                }
            },
            [&](const LinkerProfile& profile) {
                if(rawData.updated & raw_data_cat::beadPosition) {
                    meshData = createLinkerMeshData(
                        select(rawData.frameData.linkers, profile.selector, rawData.displayTypeMap),
                        profile.displaySettings
                    );
                }
            },
            [&](const AuxLineProfile& profile) {
                if(rawData.updated & raw_data_cat::compartment) {
                    meshData = createAuxLineMeshData(
                        rawData.frameData,
                        profile.displaySettings
                    );
                }
            }
        };
    };

    for(auto& eachProfileData : profileData) {
        std::visit(updateMeshData(eachProfileData.data), eachProfileData.profile);
    }

    // reset updated states
    rawData.updated = 0;
}


} // namespace medyan::visual

#endif
