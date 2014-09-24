//
//  MController.cpp
//  Cyto
//
//  Created by James Komianos on 8/4/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "MController.h"

void MController::initializeFF (MechanicsFFType& forceFields) {
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

void MController::initializeMinAlgorithms (MechanicsAlgorithm& Minimizers) {
    if (Minimizers.ConjugateGradient == "FLETCHERRIEVES") {_minimizerAlgorithms.push_back(new ConjugateGradient<FletcherRieves>() );}
    else if (Minimizers.ConjugateGradient == "POLAKRIBIERE") {_minimizerAlgorithms.push_back(new ConjugateGradient<PolakRibiere>() );}
    else {
        std::cout << "Conjugate gradient method not yet implemented. Exiting" << std::endl;
        exit(EXIT_FAILURE);
    }
}

void MController::updateBeads() {
    ///Update bead-boundary interactions (VERY INEFFICIENT)
    for(auto b : *BeadDB::Instance(BeadDBKey())) b->updateBoundaryElements();
}

void MController::updateCylinders() {
    
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