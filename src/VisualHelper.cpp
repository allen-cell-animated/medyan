#include "VisualHelper.hpp"

#include <array>
#include <functional> // cref, reference_wrapper
#include <mutex>
#include <numeric> // iota
#include <thread>
#include <utility> // move
#include <vector>

#include "MathFunctions.h"
#include "Structure/Bead.h"
#include "Structure/BranchingPoint.h"
#include "Structure/Cylinder.h"
#include "Structure/Filament.h"
#include "Structure/Linker.h"
#include "Structure/MotorGhost.h"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "SysParams.h"
#include "Visual/MeshData.hpp"
#include "Visual/Render/PathExtrude.hpp"
#include "Visual/Render/Sphere.hpp"
#include "Visual/VisualElement.hpp"
#include "Visual/Window.hpp"

namespace medyan::visual {

std::weak_ptr< VisualDisplay > vdWeak;

namespace {


void helper() {
    std::lock_guard< std::mutex > sdfvGuard(sdfv.me);

    if(auto vd = vdWeak.lock()) {
        for(auto& vp : vd->vps) {
            std::lock_guard< std::mutex > veGuard(vp.veMutex);

        }
    }
} // void helper(...)

} // namespace (anonymous)

void copySystemDataAndRunHelper(sys_data_update::FlagType update) {
    if(visual::copySystemData(sdfv, update)) {

    // Launch helper thread (may use thread pool)
    std::thread(helper).detach();
    }
}

} // namespace medyan::visual
