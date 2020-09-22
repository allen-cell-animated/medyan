#ifndef MEDYAN_Visual_FrameData_hpp
#define MEDYAN_Visual_FrameData_hpp

#include <array>
#include <cstddef>
#include <fstream>
#include <vector>

#include "OutputStruct.hpp"
#include "Util/Math/Vec.hpp"

namespace medyan::visual {


//-------------------------------------
// Data structures
//-------------------------------------

struct MembraneFrame {
    // meta data
    int id = 0;

    // triangular mesh data
    std::vector< std::array< int, 3 >> triangles;
    std::vector< mathfunc::Vec3f > vertexCoords;
};


struct DisplayFrame {
    // meta data

    // frame data
    std::vector< MembraneFrame > membranes;
};


//-------------------------------------
// Functions
//-------------------------------------

inline DisplayFrame readOneFrameDataFromOutput(
    int frameNumber,
    std::istream& is,
    std::istringstream& iss
) {
    using namespace std;

    DisplayFrame res;

    OutputStructSnapshot snapshot(frameNumber);
    snapshot.getFromOutput(is, iss);

    // membrane
    {
        const int numMembranes = snapshot.membraneStruct.size();
        res.membranes.reserve(numMembranes);

        for(const auto& m : snapshot.membraneStruct) {
            MembraneFrame mf;

            // get coordinate list
            mf.vertexCoords.reserve(m.getNumVertices());
            for(const auto& coord : m.getMembraneInfo().attributeInitializerInfo.vertexCoordinateList) {
                mf.vertexCoords.push_back(mathfunc::Vec3f(coord));
            }

            // get triangle list
            mf.triangles.reserve(m.getMembraneInfo().triangleVertexIndexList.size());
            for(const auto& t : m.getMembraneInfo().triangleVertexIndexList) {
                mf.triangles.push_back({ (int)t[0], (int)t[1], (int)t[2] });
            }

            res.membranes.push_back(std::move(mf));
        }
    }

    return res;
}

inline std::vector<DisplayFrame> readAllFrameDataFromOutput() {
    using namespace std;

    vector< DisplayFrame > res;

    {
        // File location
        string snapshotFilepath("./snapshot.traj");
        // Read snapshot
        ifstream is(snapshotFilepath);

        int curFrame = 0;

        LOG(STEP) << "Start reading " << snapshotFilepath;

        string line;
        while(true) {
            getline(is, line);
            if(!is) break;
            if(line.empty()) continue;

            ++curFrame;
            if (curFrame % 20 == 0) LOG(INFO) << "Frame " << curFrame;

            istringstream iss(line);
            res.push_back(readOneFrameDataFromOutput(curFrame, is, iss));
        }
        LOG(STEP) << "Reading complete. " << res.size() << " frames loaded.";

    }


    return res;
}


} // namespace medyan::visual

#endif
