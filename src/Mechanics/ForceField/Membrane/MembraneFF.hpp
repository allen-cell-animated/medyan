#ifndef MEDYAN_Mechanics_ForceField_Membrane_MembraneFF_hpp
#define MEDYAN_Mechanics_ForceField_Membrane_MembraneFF_hpp

#include <memory> // unique_ptr
#include <stdexcept>
#include <string>
#include <vector>

#include "Mechanics/ForceField/ForceField.h"
#include "Mechanics/ForceField/Membrane/MembraneBending.hpp"
#include "Mechanics/ForceField/Membrane/MembraneBendingHelfrich.hpp"
#include "Mechanics/ForceField/Membrane/MembraneStretchingLocal.hpp"
#include "Mechanics/ForceField/Membrane/MembraneStretchingImpl.hpp"
#include "Mechanics/ForceField/Membrane/MembraneTension.hpp"
#include "Mechanics/ForceField/Membrane/MembraneTriangleProtect.hpp"
#include "Mechanics/ForceField/Types.hpp"
#include "Util/Io/Log.hpp"

struct MembraneFFFactory {

    using CurvRequirement = ForceFieldTypes::GeometryCurvRequirement;

    struct Result {
        std::vector< std::unique_ptr< ForceField > > forceFields;
        CurvRequirement curvReq = CurvRequirement::curv;
    };

    auto operator()(
        const std::string& stretchingType,
        const std::string& tensionType,
        const std::string& bendingType
    ) const {
        using namespace std;

        Result res;

        if (stretchingType == "LOCAL_HARMONIC") {
            // In material coordinates, it is mandatory, and is applicable to
            //   non-reservior-touching membrane.
            // In normal coordinates, it cannot be modeled without more degrees of freedom.
            res.forceFields.push_back(
                std::make_unique< MembraneStretchingLocal >()
            );
        }
        else if(stretchingType == "GLOBAL_HARMONIC") {
            // Assumption: surface tension is uniform on the membrane
            // Only applicable to non-reservior-touching membrane in normal coordinates.
            LOG(WARNING) << "Currently the force field is not implemented";
            // res.forceFields.push_back(
            //     std::make_unique< MembraneStretchingGlobal >()
            // );
        }
        else if(stretchingType == "") {}
        else {
            LOG(ERROR) << "Membrane stretching FF type " << stretchingType << " is not recognized.";
            throw std::runtime_error("Membrane stretching FF type not recognized");
        }

        if(tensionType == "CONSTANT") {
            // In material coordinates, it is applicable to reservior-touching border triangles.
            // In normal corodinates, it is applicable to the whole reservior-touching membrane,
            //   assuming that surface tension is constant globally.
            LOG(WARNING) << "Currently the force field is not implemented";
            // res.forceFields.push_back(
            //     std::make_unique< MembraneTension >()
            // );
        }
        else if(tensionType == "") {}
        else {
            LOG(ERROR) << "Membrane tension FF type " << tensionType << " is not recognized.";
            throw std::runtime_error("Membrane tension FF type not recognized");
        }
        
        if (bendingType == "HELFRICH") {
            res.forceFields.emplace_back(
                new MembraneBending<MembraneBendingHelfrich>()
            );
            res.curvReq = CurvRequirement::curv;
        }
        else if(bendingType == "HELFRICH_QUADRATIC") {
            res.forceFields.emplace_back(
                new MembraneBending<MembraneBendingHelfrichQuadratic>()
            );
            res.curvReq = CurvRequirement::curv2;
        }
        else if(bendingType == "") {}
        else {
            cout << "Membrane bending FF type " << bendingType << " is not recognized." << endl;
            throw std::runtime_error("Membrane bending FF type not recognized");
        }

        // Currently, do not add the protective force field
        // res.forceFields.emplace_back(
        //     new MembraneTriangleProtect< MembraneTriangleProtectFene, true >()
        // );


        return res;
    }

};

#endif
