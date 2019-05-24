#ifndef MEDYAN_Visual_SharedData_Hpp
#define MEDYAN_Visual_SharedData_Hpp

#include <mutex>
#include <vector>

namespace visual {

namespace shared {

extern std::mutex dataMutex;

extern bool coordChanged;
extern bool indexChanged;

extern std::vector< float > vertexCoords;
extern std::vector< unsigned int > triangleVertexIndices;

}
}

#endif
