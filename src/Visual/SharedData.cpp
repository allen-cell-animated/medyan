#include "Visual/SharedData.hpp"

#ifdef VISUAL

namespace visual {
namespace shared {

std::vector< float > vertexCoords;
bool coordChanged = true;

std::vector< unsigned int > vertexIndices;
bool indexChanged = true;


bool forceChanged = true;
bool forceIndexChanged = true;
std::vector< float > arrowVertexCoords;
std::vector< unsigned int > lineVertexIndices;

} // namespace shared
} // namespace visual

#endif // VISUAL
