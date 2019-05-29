#include "Visual/SharedData.hpp"

namespace visual {
namespace shared {

std::mutex dataMutex;

bool coordChanged = true;
bool indexChanged = true;

std::vector< float > vertexCoords;
std::vector< unsigned int > triangleVertexIndices;

bool forceChanged = true;
bool forceIndexChanged = true;
std::vector< float > arrowVertexCoords;
std::vector< unsigned int > lineVertexIndices;

} // namespace shared
} // namespace visual
