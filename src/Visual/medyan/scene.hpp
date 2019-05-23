#ifndef MEDYAN_VISUAL_MEDYAN_SCENE_HPP
#define MEDYAN_VISUAL_MEDYAN_SCENE_HPP

#include "visual/common.hpp"

namespace visual {

template< bool enabled > class SceneImpl;

template<> class SceneImpl< true > {
    void printSystem() const;
};

template<> class SceneImpl< false > {
    void printSystem() const {}
};

using Scene = SceneImpl< enabled >;

} // namespace visual

#endif
