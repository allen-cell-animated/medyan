#ifndef MEDYAN_Visual_FrameData_hpp
#define MEDYAN_Visual_FrameData_hpp

#include <array>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <unordered_map>
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

struct FilamentFrame {
    // meta data
    int id = 0;
    int type = 0;

    // line segment data
    std::vector< mathfunc::Vec3f > coords;
};

// Used for line-segment-like elements, such as static linkers and motors.
struct LinkerFrame {
    // meta data
    int id = 0;
    int type = 0;

    // line data
    std::array< mathfunc::Vec3f, 2 > coords;
};

// Used for point-like elements, such as bubbles and branching points (as of Oct 3 2020).
struct PointFrame {
    // meta data
    int id = 0;
    int type = 0;

    // point data
    mathfunc::Vec3f coord {};
};


struct DisplayFrame {
    // meta data
    int serial = 0;

    // frame data
    std::vector< MembraneFrame > membranes;
    std::vector< FilamentFrame > filaments;
    std::vector< LinkerFrame >   linkers;
};

struct DisplayTypeMap {
    // linker type map
    std::vector< std::string > linkerTypeName;
    std::unordered_map< std::string, int > linkerTypeMap;
};

struct DisplayData {
    // meta data
    DisplayTypeMap displayTypeMap;

    // all frames
    std::vector< DisplayFrame > frames;
};


struct DisplayTrajectoryFileSettings {
    std::filesystem::path trajSnapshot = "./snapshot.traj";
};

//-------------------------------------
// Functions
//-------------------------------------

// Get the displayed linker type from type name.
// If the name does not exist in the map, add a new one.
inline int displayLinkerType(DisplayTypeMap& displayTypeMap, const std::string& name) {
    if(auto it = displayTypeMap.linkerTypeMap.find(name); it == displayTypeMap.linkerTypeMap.end()) {
        // The name was not found
        const int res = displayTypeMap.linkerTypeName.size();
        displayTypeMap.linkerTypeMap.insert({ name, res });
        displayTypeMap.linkerTypeName.push_back(name);
        return res;
    }
    else {
        // The name was found
        return it->second;
    }
}

inline int numFrames(const DisplayData& displayData) { return displayData.frames.size(); }

// Given an open file stream, read one frame of data for display.
//
// Inputs:
//   - frameNumber: the index of the frame being read
//   - is: the file stream, starting next line
//   - iss: the string stream of the current line
//   - displayTypeMap: the id-name maps of various types
inline DisplayFrame readOneFrameDataFromOutput(
    int                 frameNumber,
    std::istream&       is,
    std::istringstream& iss,
    DisplayTypeMap&     displayTypeMap
) {
    using namespace std;

    DisplayFrame res;

    // frame number
    res.serial = frameNumber;

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

            res.membranes.push_back(move(mf));
        }
    }

    // filament
    {
        const int numFilaments = snapshot.filamentStruct.size();
        res.filaments.reserve(numFilaments);

        for(const auto& f : snapshot.filamentStruct) {
            FilamentFrame ff;

            // get id and type
            ff.id = f.getId();
            // TODO get type

            // get coordinates
            ff.coords.reserve(f.getNumBeads());
            for(const auto& coord : f.getCoords()) {
                ff.coords.push_back(mathfunc::Vec3f(coord));
            }

            res.filaments.push_back(move(ff));
        }
    }

    // linkers
    {
        const int numLinkers =
            snapshot.linkerStruct.size() + snapshot.motorStruct.size();
        res.linkers.reserve(numLinkers);

        // To accommodate for linker/motor
        const auto addLinkerFrame = [&](const auto& linkerStruct, const string& supertypeName) {
            for(const auto& l : linkerStruct) {
                LinkerFrame lf;

                lf.type = displayLinkerType(displayTypeMap, supertypeName);
                for(int i = 0; i < l.getCoords().size(); ++i) {
                    lf.coords[i] = l.getCoords()[i];
                }

                res.linkers.push_back(move(lf));
            }
        };
        addLinkerFrame(snapshot.linkerStruct, "linker");
        addLinkerFrame(snapshot.motorStruct,  "motor");
    }

    return res;
}

inline DisplayData readAllFrameDataFromOutput(
    const DisplayTrajectoryFileSettings& inputs
) {
    using namespace std;

    DisplayData res;

    {
        // Read snapshot
        ifstream is(inputs.trajSnapshot);

        int curFrame = 0;

        LOG(STEP) << "Start reading " << inputs.trajSnapshot;

        string line;
        while(true) {
            getline(is, line);
            if(!is) break;
            if(line.empty()) continue;

            ++curFrame;
            if (curFrame % 20 == 0) LOG(INFO) << "Frame " << curFrame;

            istringstream iss(line);
            res.frames.push_back(readOneFrameDataFromOutput(curFrame, is, iss, res.displayTypeMap));
        }
        LOG(STEP) << "Reading complete. " << numFrames(res) << " frames loaded.";

    }


    return res;
}


} // namespace medyan::visual

#endif
