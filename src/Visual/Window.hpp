#ifndef MEDYAN_Visual_Window_Hpp
#define MEDYAN_Visual_Window_Hpp

#include <algorithm>
#include <array>
#include <iostream> // cout, endl
#include <stdexcept> // runtime_error
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_WRITE_STATIC
#include "stb_image_write.h"

#include "Util/Environment.hpp"
#include "Util/Io/Log.hpp"
#include "Visual/Common.hpp"
#include "Visual/Control.hpp"
#include "Visual/Gui.hpp"
#include "Visual/Playback.hpp"
#include "Visual/Shader.hpp"
#include "Visual/ShaderSrc.hpp"
#include "Visual/VisualElement.hpp"
#include "VisualSystemRawData.hpp"


// For best portability, the window signal handling could only be done from the
// main thread (due to MacOS Cocoa framework).

namespace medyan::visual {


inline void glfwError(int id, const char* description) {
    LOG(ERROR) << description;
    throw std::runtime_error("Error in GLFW environment");
}



// The RAII object for managing visualization context and window
//
// Note:
//   - Only one object of this type can be created at a time.
//   - The window signal handling should only be done from the main thread (due
//     to MacOS Cocoa framework).
class VisualContext {
public:
    struct OffscreenBufferGuard {
        bool bufferBuilt = false;
        unsigned offscreenFbo = 0;
        unsigned offscreenRbo = 0;
        unsigned offscreenColorRbo = 0;

        // destructor will deallocate the frame buffers
        ~OffscreenBufferGuard() {
            unbuild();
        }

        // build a new buffer
        void build(int width, int height) {
            if(bufferBuilt) unbuild();

            bufferBuilt = true;

            glGenFramebuffers(1, &offscreenFbo);
            glBindFramebuffer(GL_FRAMEBUFFER, offscreenFbo);

            glGenRenderbuffers(1, &offscreenRbo);
            glBindRenderbuffer(GL_RENDERBUFFER, offscreenRbo);
            glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, width, height);

            glGenRenderbuffers(1, &offscreenColorRbo);
            glBindRenderbuffer(GL_RENDERBUFFER, offscreenColorRbo);
            glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA8, width, height);

            glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, offscreenRbo);
            glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, offscreenColorRbo);

            if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
                LOG(ERROR) << "Framebuffer is not complete.";
            }
        }

        void unbuild() {
            if(bufferBuilt) {
                glDeleteFramebuffers(1, &offscreenFbo);
                glDeleteRenderbuffers(1, &offscreenRbo);
                glDeleteRenderbuffers(1, &offscreenColorRbo);
                bufferBuilt = false;
            }
        }
    };



    // Display settings and states
    DisplaySettings displaySettings;
    DisplayStates   displayStates;


    VisualContext() {
        // GLFW initializing
        LOG(DEBUG) << "Initializing GLFW";
        glfwSetErrorCallback(&glfwError);
        glfwInit();
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
        // The following line is generally useless, but it is fun to see a transparent window.
        // glfwWindowHint(GLFW_TRANSPARENT_FRAMEBUFFER, GL_TRUE);

        // Window initializing
        window_ = glfwCreateWindow(displaySettings.mainView.canvas.width, displaySettings.mainView.canvas.height, "MEDYAN", NULL, NULL);
        if(window_ == NULL) {
            LOG(ERROR) << "Failed to create GLFW window";
            glfwTerminate();
            return;
        }
        glfwMakeContextCurrent(window_);

        // GLAD initializing
        LOG(DEBUG) << "initializing GLAD";
        if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
            LOG(ERROR) << "Failed to initialize GLAD";
            return;
        }

        glViewport(0, 0, displaySettings.mainView.canvas.width, displaySettings.mainView.canvas.height);
        glfwSetWindowUserPointer(window_, this);

        // Set window callbacks
        const auto framebufferSizeCallback = [](GLFWwindow* window, int width, int height) {
            auto& vc = *static_cast< VisualContext* >(glfwGetWindowUserPointer(window));
            auto& canvas = vc.displaySettings.mainView.canvas;
            canvas.width = width;
            canvas.height = height;
            glViewport(0, 0, width, height);
        };
        const auto cursorPositionCallback = [](GLFWwindow* window, double xpos, double ypos) {
            auto& vc = *static_cast< VisualContext* >(glfwGetWindowUserPointer(window));
            auto& controlStates = vc.displayStates.mainView.control;
            auto& camera = vc.displaySettings.mainView.camera;
            const auto& control = vc.displaySettings.mainView.control;

            const int mouseState = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
            if(mouseState == GLFW_PRESS) {
                if(controlStates.mouseLeftAlreadyPressed) {
                    mouseDragCamera(
                        vc.displaySettings.mainView.camera,
                        vc.displayStates.mainView.control.mouseLastX,
                        vc.displayStates.mainView.control.mouseLastY,
                        xpos,
                        ypos,
                        vc.displaySettings.mainView.control
                    );

                } else {
                    controlStates.mouseLeftAlreadyPressed = true;
                }
                controlStates.mouseLastX = xpos;
                controlStates.mouseLastY = ypos;
            } else {
                controlStates.mouseLeftAlreadyPressed = false;
            }
        };
        const auto scrollCallback = [](GLFWwindow* window, double xoffset, double yoffset) {
            auto& vc = *static_cast< VisualContext* >(glfwGetWindowUserPointer(window));
            auto& proj = vc.displaySettings.mainView.projection;
            auto& fov = proj.fov;
            auto& scale = proj.scale;

            fov -= 0.02 * yoffset;
            fov = std::clamp(fov, 0.01f, 3.13f);

            scale -= 0.05 * yoffset;
            scale = std::clamp(scale, 0.05f, 20.0f);
        };
        const auto keyCallback = [](GLFWwindow* window, int key, int scancode, int action, int mods) {
            auto& vc = *static_cast< VisualContext* >(glfwGetWindowUserPointer(window));

            // Press F to (pay respect) take snapshot.
            if(key == GLFW_KEY_F && action == GLFW_PRESS) {
                vc.displayStates.mainView.control.snapshotRenderingNextFrame = true;
            }

            // Press P to generate snapshots for all frames.
            if(key == GLFW_KEY_P && action == GLFW_PRESS) {
                if(
                    auto& offscreenRender = vc.displayStates.playback.offscreenRender;
                    vc.displaySettings.displayMode == DisplayMode::trajectory &&
                    !offscreenRender.has_value()
                ) {
                    offscreenRender.emplace(
                        TrajectoryPlaybackStates::OffscreenRender {
                            0,
                            vc.displayStates.playback.maxFrame
                        }
                    );
                    // Reset the playhead to be the previous frame of the minimum of the range of frames,
                    // even if the frame does not exist.
                    vc.displayStates.playback.currentFrame = offscreenRender->frameRangeLo - 1;
                }
            }

            // Press G to toggle GUI.
            if(key == GLFW_KEY_G && action == GLFW_PRESS) {
                vc.displaySettings.gui.enabled = !vc.displaySettings.gui.enabled;
            }
        };

        glfwSetFramebufferSizeCallback(window_, framebufferSizeCallback);
        glfwSetCursorPosCallback(window_, cursorPositionCallback);
        glfwSetScrollCallback(window_, scrollCallback);
        glfwSetKeyCallback(window_, keyCallback);

    } // VisualContext()

    ~VisualContext() {
        glfwTerminate();
    }

    auto window() const { return window_; }

    // Helper function to process window inputs
    void processInput() {

        auto& camera = displaySettings.mainView.camera;

        const float currentTime = glfwGetTime();
        const float deltaTime = currentTime - displayStates.timing.glfwTimeLastFrame;
        displayStates.timing.update(currentTime);
        const float cameraMove = displaySettings.mainView.control.cameraKeyPositionPerFrame * deltaTime;

        if (glfwGetKey(window_, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
            // Do nothing
        }

        if(glfwGetKey(window_, GLFW_KEY_W) == GLFW_PRESS) {
            const auto change = cameraMove * glm::normalize(camera.target - camera.position);
            camera.position += change;
            camera.target += change;
        }
        if(glfwGetKey(window_, GLFW_KEY_S) == GLFW_PRESS) {
            const auto change = cameraMove * glm::normalize(camera.target - camera.position);
            camera.position -= change;
            camera.target -= change;
        }
        if(glfwGetKey(window_, GLFW_KEY_A) == GLFW_PRESS) {
            const auto change = glm::normalize(glm::cross(glm::normalize(camera.target - camera.position), camera.up)) * cameraMove;
            camera.position -= change;
            camera.target -= change;
        }
        if(glfwGetKey(window_, GLFW_KEY_D) == GLFW_PRESS) {
            const auto change = glm::normalize(glm::cross(glm::normalize(camera.target - camera.position), camera.up)) * cameraMove;
            camera.position += change;
            camera.target += change;
        }
    }

private:

    // Member variables
    GLFWwindow* window_;

}; // VisualContext


// The RAII object for all the rendering process
struct VisualDisplay {
    // The overall opengl context. Must be at top
    VisualContext vc;

    // The ImGui context. Must be below the opengl context
    ImguiGuard guard;

    // Shader objects
    Shader shaderSurface { shader::VertexElementLight, shader::FragElementLight };
    Shader shaderLine    { shader::VertexElementLine,  shader::FragElementLine  };


    VisualDisplay(DisplayMode displayMode = DisplayMode::realtime) :
        guard(vc.window())
    {
        // Configure global opengl state
        glEnable(GL_DEPTH_TEST);

        // Set initial display mode
        vc.displaySettings.displayMode = displayMode;

        // Setup realtime display profiles
        vc.displayStates.realtimeDataStates.profileData = makeDefaultElementProfileData();

    } // VisualDisplay()

    void run() {

        while (!glfwWindowShouldClose(vc.window())) {
            // input
            vc.processInput();

            // check for update for mesh data
            if(vc.displaySettings.displayMode == DisplayMode::trajectory) {
                appendTrajectoryDataFromLoadDock(vc.displayStates);
                playbackUpdateMaxFrames(vc.displayStates);

                // check whether updates are needed by new frame.
                playbackCheckTrajectory(vc.displaySettings, vc.displayStates, glfwGetTime());

                // update mesh data if needed.
                updateMeshDataForAllTrajectories(vc.displayStates.trajectoryDataStates, vc.displayStates.playback.currentFrame);
            } else {
                // update mesh data if needed.
                updateRealtimeMeshData(vc.displayStates.realtimeDataStates.profileData, sdfv);
            }

            auto& mainViewSettings = vc.displaySettings.mainView;
            auto& mainViewStates   = vc.displayStates.mainView;

            // select frame buffer: on screen display / offscreen snapshot
            const bool offscreen =
                mainViewStates.control.snapshotRenderingNextFrame ||
                vc.displayStates.playback.offscreenRender.has_value();
            mainViewStates.control.snapshotRenderingNextFrame = false;

            auto width  = mainViewSettings.canvas.width;
            auto height = mainViewSettings.canvas.height;
            auto projectionWidth  = mainViewSettings.canvas.width;
            auto projectionHeight = mainViewSettings.canvas.height;
            if(offscreen) {
                std::tie(width, height) = mainViewSettings.control.snapshotSize(width, height);
                if(
                    mainViewSettings.projection.type == ObjectViewSettings::Projection::Type::orthographic &&
                    mainViewSettings.control.snapshotResolution == ObjectViewSettings::Control::SnapshotResolution::scaleWithScreen &&
                    mainViewSettings.control.snapshotUndoScaleOrtho
                ) {
                    // The projection size stay as is.
                    // so do nothing
                } else {
                    projectionWidth = width;
                    projectionHeight = height;
                }
            }

            VisualContext::OffscreenBufferGuard offscreenBufferGuard;
            if(offscreen) {
                offscreenBufferGuard.build(width, height);
                LOG(STEP) << "Rendering offscreen...";
            }
            glBindFramebuffer(GL_FRAMEBUFFER, offscreen ? offscreenBufferGuard.offscreenFbo : 0);

            glViewport(0, 0, width, height);
            {
                const auto& c = mainViewSettings.canvas.bgColor;
                glClearColor(c[0], c[1], c[2], c[3]);
            }
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            // transform
            const auto model          = glm::mat4(1.0f);
            const auto modelInvTrans3 = glm::transpose(glm::inverse(model));
            const auto view           = mainViewSettings.camera.view();
            const auto projection     = mainViewSettings.projection.proj(projectionWidth, projectionHeight);

            // Start to render surface
            glUseProgram(shaderSurface.id());

            shaderSurface.setMat4("projection",     projection);
            shaderSurface.setMat4("model",          model);
            shaderSurface.setMat3("modelInvTrans3", modelInvTrans3);
            shaderSurface.setMat4("view",           view);

            shaderSurface.setVec3("CameraPos",   mainViewSettings.camera.position);

            shaderSurface.setVec3("dirLights[0].direction", glm::vec3 {1.0f, 1.0f, 1.0f});
            shaderSurface.setVec3("dirLights[0].ambient",   glm::vec3 {0.1f, 0.1f, 0.1f});
            shaderSurface.setVec3("dirLights[0].diffuse",   glm::vec3 {0.3f, 0.3f, 0.3f});
            shaderSurface.setVec3("dirLights[0].specular",  glm::vec3 {0.5f, 0.5f, 0.5f});
            shaderSurface.setVec3("dirLights[1].direction", glm::vec3 {-1.0f, -1.0f, -1.0f});
            shaderSurface.setVec3("dirLights[1].ambient",   glm::vec3 {0.1f, 0.1f, 0.1f});
            shaderSurface.setVec3("dirLights[1].diffuse",   glm::vec3 {0.3f, 0.3f, 0.3f});
            shaderSurface.setVec3("dirLights[1].specular",  glm::vec3 {0.5f, 0.5f, 0.5f});

            const glm::vec3 pointLightPositions[4] {
                { -500.0f, -500.0f, -500.0f },
                { -500.0f, 3500.0f, 3500.0f },
                { 3500.0f, -500.0f, 3500.0f },
                { 3500.0f, 3500.0f, -500.0f }
            };
            shaderSurface.setVec3("pointLights[0].position", pointLightPositions[0]);
            shaderSurface.setVec3("pointLights[0].ambient",  glm::vec3 { 0.05f, 0.05f, 0.05f });
            shaderSurface.setVec3("pointLights[0].diffuse",  glm::vec3 { 0.6f, 0.6f, 0.6f });
            shaderSurface.setVec3("pointLights[0].specular", glm::vec3 { 1.0f, 1.0f, 1.0f });
            shaderSurface.setFloat("pointLights[0].constant",  1.0f);
            shaderSurface.setFloat("pointLights[0].linear",    1.4e-4f);
            shaderSurface.setFloat("pointLights[0].quadratic", 7.2e-8f);
            shaderSurface.setVec3("pointLights[1].position", pointLightPositions[1]);
            shaderSurface.setVec3("pointLights[1].ambient",  glm::vec3 { 0.05f, 0.05f, 0.05f });
            shaderSurface.setVec3("pointLights[1].diffuse",  glm::vec3 { 0.6f, 0.6f, 0.6f });
            shaderSurface.setVec3("pointLights[1].specular", glm::vec3 { 1.0f, 1.0f, 1.0f });
            shaderSurface.setFloat("pointLights[1].constant",  1.0f);
            shaderSurface.setFloat("pointLights[1].linear",    1.4e-4f);
            shaderSurface.setFloat("pointLights[1].quadratic", 7.2e-8f);
            shaderSurface.setVec3("pointLights[2].position", pointLightPositions[2]);
            shaderSurface.setVec3("pointLights[2].ambient",  glm::vec3 { 0.05f, 0.05f, 0.05f });
            shaderSurface.setVec3("pointLights[2].diffuse",  glm::vec3 { 0.1f, 0.1f, 0.1f });
            shaderSurface.setVec3("pointLights[2].specular", glm::vec3 { 0.2f, 0.2f, 0.2f });
            shaderSurface.setFloat("pointLights[2].constant",  1.0f);
            shaderSurface.setFloat("pointLights[2].linear",    1.4e-4f);
            shaderSurface.setFloat("pointLights[2].quadratic", 7.2e-8f);
            shaderSurface.setVec3("pointLights[3].position", pointLightPositions[3]);
            shaderSurface.setVec3("pointLights[3].ambient",  glm::vec3 { 0.05f, 0.05f, 0.05f });
            shaderSurface.setVec3("pointLights[3].diffuse",  glm::vec3 { 0.1f, 0.1f, 0.1f });
            shaderSurface.setVec3("pointLights[3].specular", glm::vec3 { 0.2f, 0.2f, 0.2f });
            shaderSurface.setFloat("pointLights[3].constant",  1.0f);
            shaderSurface.setFloat("pointLights[3].linear",    1.4e-4f);
            shaderSurface.setFloat("pointLights[3].quadratic", 7.2e-8f);

            const auto drawSurfaceProfileData = [&](auto& eachProfileData) {
                if(displayGeometryType(eachProfileData.profile) == DisplayGeometryType::surface) {

                    std::visit(
                        [&](const auto& profile) {
                            draw(eachProfileData.data, std::nullopt, shaderSurface, profile.displaySettings);
                        },
                        eachProfileData.profile
                    );
                }
            };
            if(vc.displaySettings.displayMode == DisplayMode::trajectory) {
                for(auto& traj : vc.displayStates.trajectoryDataStates.trajectories) {
                    if(traj.displayMasterSwitch) {
                        for(auto& profileData : traj.profileData) {
                            drawSurfaceProfileData(profileData);
                        }
                    }
                }
            }
            else {
                for(auto& profileData : vc.displayStates.realtimeDataStates.profileData) {
                    drawSurfaceProfileData(profileData);
                }
            }


            // start to render line
            glUseProgram(shaderLine.id());

            shaderLine.setMat4("projection",     projection);
            shaderLine.setMat4("model",          model);
            shaderLine.setMat3("modelInvTrans3", modelInvTrans3);
            shaderLine.setMat4("view",           view);

            const auto drawLineProfileData = [&](auto& eachProfileData) {
                if(displayGeometryType(eachProfileData.profile) == DisplayGeometryType::line) {

                    std::visit(
                        [&](const auto& profile) {
                            draw(eachProfileData.data, std::nullopt, shaderLine, profile.displaySettings);
                        },
                        eachProfileData.profile
                    );
                }
            };
            if(vc.displaySettings.displayMode == DisplayMode::trajectory) {
                for(auto& traj : vc.displayStates.trajectoryDataStates.trajectories) {
                    for(auto& profileData : traj.profileData) {
                        drawLineProfileData(profileData);
                    }
                }
            }
            else {
                for(auto& profileData : vc.displayStates.realtimeDataStates.profileData) {
                    drawLineProfileData(profileData);
                }
            }


            glBindVertexArray(0);

            // Output offscreen render results
            if(offscreen) {
                std::vector< std::uint8_t > data(width * height * 4); // 8-bit color
                glReadBuffer(GL_COLOR_ATTACHMENT0);
                glReadPixels(0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, data.data());

                int frameIndex = -1;
                if(vc.displayStates.playback.offscreenRender.has_value()) {
                    frameIndex = vc.displayStates.playback.currentFrame;
                }
                const auto snapshotPngFile = mainViewSettings.control.snapshotPngFile(frameIndex);
                stbi_write_png(snapshotPngFile.string().c_str(), width, height, 4, data.data(), 4 * width);
                LOG(INFO) << "Snapshot saved to " << snapshotPngFile;
            }

            // Update GUI
            imguiLoopRender(vc.displaySettings, vc.displayStates);
            // Perform async tasks
            dispatchAnAsyncTask(vc.displayStates.sync);

            // check
            glfwSwapBuffers(vc.window());
            glfwPollEvents();

        } // End main loop

    } // void run() const
}; // VisualDisplay

} // namespace medyan::visual

#endif
