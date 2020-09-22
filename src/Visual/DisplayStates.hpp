#ifndef MEDYAN_Visual_DisplayStates_hpp
#define MEDYAN_Visual_DisplayStates_hpp

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
};

} // namespace medyan::visual

#endif
