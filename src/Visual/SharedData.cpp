#include "Visual/SharedData.hpp"

namespace visual {
namespace shared {

std::mutex dataMutex;

bool coordChanged = true;
bool indexChanged = true;

std::vector< float > vertexCoords;
std::vector< unsigned int > triangleVertexIndices;

} // namespace shared
} // namespace visual
