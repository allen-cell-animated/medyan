#include "Visual/SharedData.hpp"

namespace visual {
namespace shared {

std::mutex dataMutex;

std::vector< std::shared_ptr< VisualElement > > visualElements;

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
