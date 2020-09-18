#ifndef MEDYAN_Visual_DisplaySettings_hpp
#define MEDYAN_Visual_DisplaySettings_hpp

namespace medyan::visual {

enum class DisplayMode {
    realtime,        // Display system data from memory during simulation
    trajectory       // Read and display trajectory from output files
};

struct DisplaySettings {
    DisplayMode displayMode = DisplayMode::realtime;

};

} // namespace medyan::visual

#endif
