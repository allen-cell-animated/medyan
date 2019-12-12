#ifndef MEDYAN_VisualHelper_Hpp
#define MEDYAN_VisualHelper_Hpp

#include <cstdint>
#include <memory> // weak_ptr

#include "VisualSystemRawData.hpp"

//-----------------------------------------------------------------------------
// Visualization helper transforms the system related data to actual OpenGL
// recognized data according to the settings.
//
// Note: the helper function should only be called from the simulation thread.
// This could be either called directly (serial) or launched (and detached) in
// a different thread
//-----------------------------------------------------------------------------

namespace visual {

struct VisualDisplay;
extern std::weak_ptr< VisualDisplay > vdWeak;

namespace sys_data_update {

    using FlagType = std::uint_fast8_t;

    constexpr FlagType None            = 0;
    constexpr FlagType BeadPosition    = 1 << 0;
    constexpr FlagType BeadConnection  = 1 << 1;
    constexpr FlagType Compartment     = 1 << 2;
    constexpr FlagType Concentration   = 1 << 3;

} // namespace sys_data_update

void copySystemDataAndRunHelper(sys_data_update::FlagType update);

} // namespace visual

#endif
