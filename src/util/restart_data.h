#ifndef MEDYAN_UTIL_RESTART_DATA_H
#define MEDYAN_UTIL_RESTART_DATA_H

#include <array>
#include <cstdint>
#include <vector>

#include "util/types.h"

namespace restart {

// RestartData contains sufficient data for rebuilding the whole simulation,
// and can be used by the system and restart_io processes.
struct RestartData {
    // Types
    struct FilamentData {
        std::uint16_t type;
        std::vector< std::array< mreal, 3 > > beadPosition;
    };

    std::vector< FilamentData > filamentData;
};

} // namespace restart

#endif
