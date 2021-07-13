#ifndef MEDYAN_Core_Globals_Hpp
#define MEDYAN_Core_Globals_Hpp

#include <string>

namespace medyan {

struct GlobalVar {
    // Analyze specific
    std::size_t analyzeMembraneBondFrame = 0;
    bool        analyzeMembraneBondAllFrames = false; // If it is set to true, it overrides the specific bond frame info.
    int analyzeFrameInterval = 1; // must be >= 1

};

inline auto& globalMutable() { static GlobalVar g; return g; }
inline const auto& global() { return globalMutable(); }

} // namespace medyan


#endif