#ifndef MEDYAN_Visual_DisplayStates_hpp
#define MEDYAN_Visual_DisplayStates_hpp

#include <atomic>
#include <functional>
#include <queue>
#include <thread>

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
    };

    Control control;
};

struct SyncStates {
    // auxiliary struct for atomic bool guard
    struct AtomicBoolGuard {
        std::atomic_bool& whichBool;
        AtomicBoolGuard(std::atomic_bool& whichBool) : whichBool(whichBool) { whichBool.store(true); }
        ~AtomicBoolGuard() { whichBool.store(false); }
    };

    // the task queue. not thread-safe
    std::queue< std::function< void() > > tasks;

    // the boolean states set by the tasks

    // at most one task can be running at a time, which will mark the busy flag true
    std::atomic_bool busy { false };

    // the individual states for discretion of the main display thread
    std::atomic_bool trajectoryLoad { false };

};

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
};

} // namespace medyan::visual

#endif
