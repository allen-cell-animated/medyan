//
//  MController.h
//  Cyto
//
//  Created by James Komianos on 8/4/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__MController__
#define __Cyto__MController__

#include "common.h"
#include "SubSystem.h"
#include "Parser.h"

#include "ForceFieldManager.h"
#include "Minimizer.h"

#include "ForceField.h"
#include "FilamentFF.h"
#include "LinkerFF.h"
#include "BoundaryFF.h"
#include "MotorGhostFF.h"

#include "ConjugateGradient.h"

#include "FilamentDB.h"

#include <iostream>
#include <vector>

class PolakRibiere;
class FletcherRieves;


/// MController class is used to initialize and run the mechanical components of a simulation

/*!
 *  MechanicalController is a class used by the SubSystem to initialize force fields, given an initial
 *  selection of which force fields should be included. It can compute forces and energies
 *  over all force fields, as well as run energy minimization algorithms.
 */

class MController {
    
private:
    ForceFieldManager _FFManager;  ///<container and methods for all force fields in system
    std::vector<Minimizer*> _minimizerAlgorithms; ///<vector with algorythms for system equlibration
    SubSystem* _subSystem;
    
    ///Initialize the MController using a list of vector names
    ///@param forceFields - a list of forcefields to be added
    void initializeFF (MechanicsFFType forceFields)
    {
        /// Check if exist!!!
        _FFManager._forceFields.push_back(new FilamentFF(forceFields.FStretchingType, forceFields.FBendingType, forceFields.FTwistingType) );
        std::cout << "Filament force field initialized: " <<std::endl;
        if(forceFields.FStretchingType != "") std::cout << "Stretching: " << forceFields.FStretchingType << std::endl;
        if(forceFields.FBendingType != "") std::cout << "Bending: " << forceFields.FBendingType<< std::endl;
        if(forceFields.FTwistingType != "")std::cout << "Twisting: " << forceFields.FTwistingType <<std::endl;
        
        _FFManager._forceFields.push_back(new LinkerFF(forceFields.LStretchingType, forceFields.LBendingType, forceFields.LTwistingType) );
        std::cout << "Linker force field initialized:"<<std::endl;
        if(forceFields.LStretchingType != "") std::cout << "Stretching: " << forceFields.LStretchingType<< std::endl;
        if(forceFields.LBendingType != "") std::cout << "Bending: " << forceFields.LBendingType<< std::endl;
        if(forceFields.LTwistingType != "") std::cout << "Twisting: " << forceFields.LTwistingType <<std::endl;
        
        _FFManager._forceFields.push_back(new MotorGhostFF(forceFields.MStretchingType, forceFields.MBendingType, forceFields.MTwistingType) );
        std::cout << "Motor force field initialized:"<<std::endl;
        if(forceFields.MStretchingType != "") std::cout << "Stretching: " << forceFields.MStretchingType<< std::endl;
        if(forceFields.MBendingType != "") std::cout << "Bending: " << forceFields.MBendingType<< std::endl;
        if(forceFields.MTwistingType != "") std::cout << "Twisting: " << forceFields.MTwistingType <<std::endl;
        
        _FFManager._forceFields.push_back(new BoundaryFF(forceFields.BoundaryFFType, "", "") );
        std::cout << "Boundary force field initialized: " <<std::endl;
        if(forceFields.BoundaryFFType != "") std::cout << "Boundary: " << forceFields.BoundaryFFType << std::endl;
        
    }
    
    void initializeMinAlgorithms (MechanicsAlgorithm Minimizers) {
        if (Minimizers.ConjugateGradient == "FLETCHERRIEVES") {_minimizerAlgorithms.push_back(new ConjugateGradient<FletcherRieves>() );}
        else if (Minimizers.ConjugateGradient == "POLAKRIBIERE") {_minimizerAlgorithms.push_back(new ConjugateGradient<PolakRibiere>() );}
        else {
            std::cout << "Conjugate gradient method not yet implemented. Exiting" << std::endl;
            exit(EXIT_FAILURE);
        }
    }

public:
    
    MController(SubSystem* s) {_subSystem = s;}
    
    void initialize(MechanicsFFType forceFields, MechanicsAlgorithm Minimizers )
    {
        initializeFF(forceFields);
        initializeMinAlgorithms(Minimizers);
    }

    
    ///Run minimization on the system using the chosen algorithm
    void run() {

        _minimizerAlgorithms[0]->Equlibrate(_FFManager);
        
        ///Update bead-boundary interactions (VERY INEFFICIENT)
        for(auto b : *BeadDB::Instance(BeadDBKey())) b->updateBoundaryElements();
        
#ifdef CHEMISTRY
        ///Update cylinder positions (ALSO VERY INEFFICIENT)
        for(auto &f : *FilamentDB::Instance(FilamentDBKey())) {
            
            for(auto &Cyl : f->getCylinderVector()) {
                std::vector<double> midpoint = mathfunc::MidPointCoordinate(Cyl->getMCylinder()->GetFirstBead()->coordinate,
                                                                            Cyl->getMCylinder()->GetSecondBead()->coordinate,
                                                                            0.5);
                Compartment* c;
                try {c = GController::getCompartment(midpoint);}
                catch (std::exception& e) {std:: cout << e.what(); exit(EXIT_FAILURE);}
                
                CCylinder* cCyl = Cyl->getCCylinder();
                if(c != cCyl->getCompartment()) {
                    CCylinder* clone = cCyl->clone(c);
                    Cyl->setCCylinder(clone);
                }
            }
        }
#endif
    }
    
};

#endif /* defined(__Cyto__MController__) */
