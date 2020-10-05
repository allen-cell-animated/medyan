#ifndef MEDYAN_Visual_DisplayStates_hpp
#define MEDYAN_Visual_DisplayStates_hpp

#include <atomic>
#include <functional>
#include <mutex>
#include <queue>
#include <thread>

#include "Visual/FrameData.hpp"
#include "Visual/VisualElement.hpp"

// Unlike DisplaySettings, the DisplayStates structure contains states for
// display thread that is generally automatically filled internally.
//
// It may also contain synchronization elements that cannot be easily copied
// or moved

namespace medyan::visual {

struct ObjectViewStates {
    struct Control {
        // mouse states
        bool mouseLeftAlreadyPressed = false;
        double mouseLastX = 0;
        double mouseLastY = 0;

        // framebuffer control
        bool snapshotRenderingNextFrame = false;
    };

    Control control;
};

struct SyncStates {

    //---------------------------------
    // Implementation of one background task
    //---------------------------------

    // the task queue. not thread-safe
    std::queue< std::function< void() > > tasks;

    // the boolean states set by the tasks

    // at most one task can be running at a time, which will mark the busy flag true
    std::atomic_bool busy { false };

    // the individual states for discretion of the main display thread
    std::mutex  meTrajectoryLoad;
    std::queue<
        std::tuple< DisplayTrajectoryFileSettings, DisplayData >
    >           trajectoryLoadDock;             // guarded by meTrajectoryLoad

};


struct DisplayTrajectoryDataStates {
    struct Trajectory {
        DisplayTrajectoryFileSettings inputs;
        DisplayData                   data;

        std::vector< ProfileWithMeshData > profileData;

        Trajectory(Trajectory&&) = default;
    };

    std::vector< Trajectory > trajectories;

};

struct DisplayRealtimeDataStates {
    std::vector< ProfileWithMeshData > profileData;
};

struct TrajectoryPlaybackStates {
    // Play back states
    int currentFrame = 0;
    int previousFrame = 0;
    int maxFrame = 0;
    bool isPlaying = true;
    float lastGlfwTime = 0;

    bool meshDataGenerated = false; // set true by mesh data generation, set false by changing frame
};


struct DisplayStates {

    struct Timing {
        // settings
        float fpsUpdateTimeInterval = 1;

        // states
        float glfwTimeLastFrame = 0;

        float glfwTimeLastFpsCheckpoint = 0;
        int   numFramesSinceLastFpsCheckpoint = 0;
        float fps = 0;

        // functions

        // Update with the time of new frame, and calculate fps
        void update(float glfwTimeNewFrame) {
            // update fps
            const float deltaCheckpoint = glfwTimeNewFrame - glfwTimeLastFpsCheckpoint;
            ++numFramesSinceLastFpsCheckpoint;
            if(deltaCheckpoint > fpsUpdateTimeInterval) {
                fps = numFramesSinceLastFpsCheckpoint / deltaCheckpoint;
                glfwTimeLastFpsCheckpoint = glfwTimeNewFrame;
                numFramesSinceLastFpsCheckpoint = 0;
            }

            // update last frame time
            glfwTimeLastFrame = glfwTimeNewFrame;
        }
    };

    // timer states
    Timing timing;

    // main view states
    ObjectViewStates mainView;

    // task and synchronization
    SyncStates sync;

    // The trajectory data states
    DisplayTrajectoryDataStates trajectoryDataStates;
    // The realtime data profiles
    DisplayRealtimeDataStates realtimeDataStates;

    // Playback states
    TrajectoryPlaybackStates playback;

};


//-------------------------------------
// Functions
//-------------------------------------

// Push a task to the background task queue, which will be serially executed.
// The task will unset the "busy" flag upon exiting.
//
// Notes:
//   - F must be callable with signature void()
template< typename F >
inline void pushAnAsyncTask(SyncStates& sync, F&& f) {
    sync.tasks.push([&sync, f] {
        // set sync.busy to false on exit
        struct BusyFalseGuard {
            std::atomic_bool& busy;
            BusyFalseGuard(std::atomic_bool& busy) : busy(busy) {}
            ~BusyFalseGuard() { busy.store(false); }
        };
        BusyFalseGuard bfg(sync.busy);

        // execute the task
        f();
    });
}

// Dispatch the first element in the background task queue.
// Before executing, set the "busy" flag.
inline void dispatchAnAsyncTask(SyncStates& sync) {
    if(!sync.tasks.empty()) {
        // if not busy, set the state to busy and execute the task
        // otherwise, do nothing
        bool expectedBusy = false;
        if(sync.busy.compare_exchange_strong(expectedBusy, true)) {

            std::thread exe(std::move(sync.tasks.front()));
            sync.tasks.pop();
            exe.detach();
        }
    }
}


// Background task: read trajectory into DisplayData
inline void backgroundTaskReadTrajectory(
    DisplayStates& states,
    const DisplayTrajectoryFileSettings& inputs
) {
    // Read trajectory
    auto displayData = readAllFrameDataFromOutput(inputs);

    // Add to loading dock
    std::scoped_lock lk(states.sync.meTrajectoryLoad);

    states.sync.trajectoryLoadDock.push(std::tuple {
        inputs,
        std::move(displayData)
    });
}

// Visual thread task: append trajectory data from loading dock
inline void appendTrajectoryDataFromLoadDock(
    DisplayStates& states
) {
    std::unique_lock lk(states.sync.meTrajectoryLoad, std::try_to_lock);
    if(!lk.owns_lock()) {
        return;
    }

    auto& dock = states.sync.trajectoryLoadDock;
    while(!dock.empty()) {
        auto [inputs, displayData] = std::move(dock.front());
        dock.pop();

        states.trajectoryDataStates.trajectories.push_back({
            std::move(inputs),
            std::move(displayData),
            makeDefaultElementProfileData()
        });
    }
}

} // namespace medyan::visual

#endif
