#ifndef MEDYAN_Util_Parser_MembraneParser_hpp
#define MEDYAN_Util_Parser_MembraneParser_hpp

#include "Parser.h"
#include "SysParams.h"

namespace medyan {

inline const auto& membraneSetupParser() {
    using namespace std;

    static auto res = [] {
        KeyValueParser< MembraneSetup > p;

        p.addStringArgs(
            "vertex-system",
            [](MembraneSetup& ms, const vector<string>& input) {
                if(input.size() != 2) {
                    LOG(ERROR) << "Only one parameter is allowed for vertex-system.";
                    throw std::runtime_error("Invalid argument.");
                }
                auto& v = input[1];
                if(v == "material") {
                    ms.vertexSystem = MembraneMeshVertexSystem::material;
                }
                else if(v == "normal") {
                    ms.vertexSystem = MembraneMeshVertexSystem::normal;
                }
                else if(v == "general") {
                    ms.vertexSystem = MembraneMeshVertexSystem::general;
                }
                else {
                    LOG(ERROR) << "Invalid vertex-system: " << v;
                    throw std::runtime_error("Invalid argument.");
                }
            },
            [](const MembraneSetup& ms) {
                vector<string> res;
                switch(ms.vertexSystem) {
                    case MembraneMeshVertexSystem::material:
                        res.push_back("material");
                        break;
                    case MembraneMeshVertexSystem::normal:
                        res.push_back("normal");
                        break;
                    case MembraneMeshVertexSystem::general:
                        res.push_back("general");
                        break;
                    default:
                        LOG(ERROR) << "Unknown membrane mesh vertex system.";
                        res.push_back("error");
                }
                return res;
            }
        );
        p.addSingleArg("area-k",    [](auto&& ms) -> auto& { return ms.areaElasticity; });
        p.addSingleArg("bending-k", [](auto&& ms) -> auto& { return ms.bendingElasticity; });
        p.addSingleArg("eq-curv",   [](auto&& ms) -> auto& { return ms.eqMeanCurv; });
        p.addSingleArg("tension",   [](auto&& ms) -> auto& { return ms.tension; });
        p.addSingleArg("volume-k",  [](auto&& ms) -> auto& { return ms.volumeElasticity; });

        p.addSingleArg("init-eq-area-factor", [](auto&& ms) -> auto& { return ms.initEqAreaFactor; });
        p.addStringArgs(
            "init",
            [](MembraneSetup& ms, const vector<string>& input) {
                auto newParam = input;
                newParam.erase(newParam.begin());
                ms.meshParam.push_back(move(newParam));
            },
            [](const MembraneSetup& ms) {
                return ms.meshParam;
            }
        );

        return p;
    }();

    return res;
}

} // namespace medyan

#endif
