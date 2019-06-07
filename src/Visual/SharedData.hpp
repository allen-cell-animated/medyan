#ifndef MEDYAN_Visual_SharedData_Hpp
#define MEDYAN_Visual_SharedData_Hpp

#include <memory> // shared_ptr
#include <mutex>
#include <vector>

#include "Visual/VisualElement.hpp"

namespace visual {
namespace shared {

extern std::mutex dataMutex;

extern std::mutex veMutex;
// Must be created and deleted within a valid opengl context
extern std::vector< std::shared_ptr< VisualElement > > visualElements;

extern std::vector< float > vertexCoords;
extern bool coordChanged;

extern std::vector< unsigned int > vertexIndices;
extern bool indexChanged;

} // namespace shared
} // namespace visual

#endif
