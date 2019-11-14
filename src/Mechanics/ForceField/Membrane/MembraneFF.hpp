#ifndef MEDYAN_Mechanics_ForceField_Membrane_MembraneFF_hpp
#define MEDYAN_Mechanics_ForceField_Membrane_MembraneFF_hpp

#include <memory> // unique_ptr
#include <stdexcept>
#include <string>
#include <vector>

#include "Mechanics/ForceField/ForceField.h"
#include "Mechanics/ForceField/Membrane/MembraneBending.hpp"
#include "Mechanics/ForceField/Membrane/MembraneBendingHelfrich.hpp"
#include "Mechanics/ForceField/Membrane/MembraneStretching.hpp"
#include "Mechanics/ForceField/Membrane/MembraneStretchingImpl.hpp"
#include "Mechanics/ForceField/Membrane/MembraneTriangleProtect.hpp"
#include "Util/Io/Log.hpp"

struct MembraneFFFactory {

    enum class CurvatureRequirement {
        curv, curv2
    };

    struct Result {
        std::vector< std::unique_ptr< ForceField > > forceFields;
        CurvatureRequirement curvReq;
    };

    auto operator()(
        const std::string& stretchingType,
        const std::string& stretchingAccuType,
        const std::string& bendingType
    ) const {
        using namespace std;

        Result res;

        if (stretchingType == "HARMONIC")
            if(stretchingAccuType == "TRIANGLE")
                res.forceFields.emplace_back(
                    new MembraneStretching< MembraneStretchingHarmonic, MembraneStretchingAccumulationType::ByTriangle >()
                );
            else if(stretchingAccuType == "VERTEX")
                res.forceFields.emplace_back(
                    new MembraneStretching< MembraneStretchingHarmonic, MembraneStretchingAccumulationType::ByVertex >()
                );
            else {
                LOG(ERROR) << "Membrane stretching accumulation type " << stretchingAccuType << " is not recognized.";
                throw std::runtime_error("Membrane stretching accumulation type not recognized");
            }

        else if(stretchingType == "LINEAR")
            if(stretchingAccuType == "TRIANGLE")
                res.forceFields.emplace_back(
                    new MembraneStretching< MembraneStretchingLinear, MembraneStretchingAccumulationType::ByTriangle >()
                );
            else if(stretchingAccuType == "VERTEX")
                res.forceFields.emplace_back(
                    new MembraneStretching< MembraneStretchingLinear, MembraneStretchingAccumulationType::ByVertex >()
                );
            else {
                LOG(ERROR) << "Membrane stretching accumulation type " << stretchingAccuType << " is not recognized.";
                throw std::runtime_error("Membrane stretching accumulation type not recognized");
            }

        else if(stretchingType == "") {}

        else {
            LOG(ERROR) << "Membrane stretching FF type " << stretchingType << " is not recognized.";
            throw std::runtime_error("Membrane stretching FF type not recognized");
        }
        
        if (bendingType == "HELFRICH")
            res.forceFields.emplace_back(
                new MembraneBending<MembraneBendingHelfrich>()
            );
        else if(bendingType == "HELFRICH_QUADRATIC")
            res.forceFields.emplace_back(
                new MembraneBending<MembraneBendingHelfrichQuadratic>()
            );
        else if(bendingType == "") {}
        else {
            cout << "Membrane bending FF type " << bendingType << " is not recognized." << endl;
            throw std::runtime_error("Membrane bending FF type not recognized");
        }

        // Always add the protective force field
        res.forceFields.emplace_back(
            new MembraneTriangleProtect< MembraneTriangleProtectFene, true >()
        );


        return res;
    }

};

#endif
