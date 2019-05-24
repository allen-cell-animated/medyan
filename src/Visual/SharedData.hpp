#ifndef MEDYAN_Visual_SharedData_Hpp
#define MEDYAN_Visual_SharedData_Hpp

#include <mutex>
#include <vector>

namespace visual {

namespace shared {

std::mutex dataMutex;

bool coordChanged = true;
bool indexChanged = true;

std::vector< float > vertexCoords;
std::vector< unsigned int > triangleVertexIndices;

}
}

#endif
