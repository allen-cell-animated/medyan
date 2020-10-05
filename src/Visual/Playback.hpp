#ifndef MEDYAN_Visual_Playback_hpp
#define MEDYAN_Visual_Playback_hpp

#include <algorithm>

#include "Visual/DisplaySettings.hpp"
#include "Visual/DisplayStates.hpp"

namespace medyan::visual {

// Update max number of frames of the playback
inline void playbackUpdateMaxFrames(
    DisplayStates& states
) {
    int maxFrames = 0;
    for(const auto& traj : states.trajectoryDataStates.trajectories) {
        maxFrames = std::max<int>(maxFrames, traj.data.frames.size());
    }

    states.playback.maxFrame = max(0, maxFrames - 1);
}


// check whether a new frame should be issued
template< typename NewFrameFunc >
inline void playbackNewFrame(
    const TrajectoryPlaybackSettings& settings,
    TrajectoryPlaybackStates&         states,
    float                             glfwTime,
    NewFrameFunc&&                    func
) {
    if(states.isPlaying) {
        if(glfwTime >= states.lastGlfwTime + 1.0f / settings.fps) {
            if(states.currentFrame < states.maxFrame) ++states.currentFrame;
        }
    }

    if(states.currentFrame != states.previousFrame) {
        states.meshDataGenerated = false;

        // Execute functions for the new frame
        func();

        states.lastGlfwTime = glfwTime;
        states.previousFrame = states.currentFrame;
    }
}

inline void playbackCheckTrajectory(
    const DisplaySettings& settings,
    DisplayStates&         states,
    float                  glfwTime
) {
    playbackNewFrame(
        settings.playback,
        states.playback,
        glfwTime,
        [&] {
            for(auto& traj : states.trajectoryDataStates.trajectories) {
                for(auto& profileData : traj.profileData) {
                    if(states.playback.currentFrame < traj.data.frames.size()) {
                        profileData.data = createMeshData(
                            traj.data.frames[states.playback.currentFrame],
                            traj.data.displayTypeMap,
                            profileData.profile
                        );
                    } else {
                        profileData.data.data.clear();
                    }
                }
            }
            states.playback.meshDataGenerated = true;
        }
    );
}

} // namespace medyan::visual

#endif
