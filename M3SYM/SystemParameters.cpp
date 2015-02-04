
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#include "SystemParameters.h"

MechanicsParameters SystemParameters::MParams;
ChemistryParameters SystemParameters::CParams;
GeometryParameters SystemParameters::GParams;
BoundaryParameters SystemParameters::BParams;
DynamicRateParameters SystemParameters::DRParams;

bool SystemParameters::checkChemParameters(ChemistryData& chem) {
    
    //Check all species for consistency
    if(chem.speciesBulk.size() != CParams.numBulkSpecies) {
        
        cout << "Number of bulk species in chemistry input does not "
             << "match the system file. Check these parameters. Exiting."
             <<endl;
        
        return false;
    }
    if(chem.speciesDiffusing.size() != CParams.numDiffusingSpecies) {
        
        cout << "Number of diffusing species in chemistry input does not "
             << "match the system file. Check these parameters. Exiting."
             <<endl;
        
        return false;
    }
    if(chem.speciesFilament.size() != CParams.numFilamentSpecies) {
        
        cout << "Number of filament species in chemistry input does not "
        << "match the system file. Check these parameters. Exiting."
        <<endl;
        
        return false;
    }
    if(chem.speciesPlusEnd.size() != CParams.numPlusEndSpecies) {
        
        cout << "Number of plus end species in chemistry input does not "
             << "match the system file. Check these parameters. Exiting."
             <<endl;
        
        return false;
    }
    if(chem.speciesMinusEnd.size() != CParams.numMinusEndSpecies) {
        
        cout << "Number of minus end species in chemistry input does not "
             << "match the system file. Check these parameters. Exiting."
             <<endl;
        
        return false;
    }
    if(chem.speciesBound.size() != CParams.numBoundSpecies) {
        
        cout << "Number of bound species in chemistry input does not "
             << "match the system file. Check these parameters. Exiting."
             <<endl;
        
        return false;
    }
    if(chem.speciesLinker.size() != CParams.numLinkerSpecies) {
        
        cout << "Number of linker species in chemistry input does not "
             << "match the system file. Check these parameters. Exiting."
             <<endl;
        
        return false;
    }
    if(chem.speciesMotor.size() != CParams.numMotorSpecies) {
        
        cout << "Number of motor species in chemistry input does not "
             << "match the system file. Check these parameters. Exiting."
             <<endl;
        
        return false;
    }
    if(chem.speciesBrancher.size() != CParams.numBrancherSpecies) {
        
        cout << "Number of brancher species in chemistry input does not "
             << "match the system file. Check these parameters. Exiting."
             <<endl;
        
        return false;
    }
    //plus and minus end consistency
    if(CParams.numPlusEndSpecies < CParams.numFilamentSpecies) {
        cout << "There must be a plus end for every filament species. Exiting."
             << endl;
        return false;
    }
    if(CParams.numMinusEndSpecies < CParams.numFilamentSpecies) {
        cout << "There must be a minus end for every filament species. Exiting."
             << endl;
        exit(EXIT_FAILURE);
    }
    
    return true;
}

bool SystemParameters::checkMechParameters(MechanicsFFType& mech) {
    
    //check ff and associated parameters for consistency
    
    //FILAMENT
    if(mech.FStretchingType != "" && MParams.FStretchingK == 0) {
        cout << "Must set a filament stetching constant. Exiting." << endl;
        return false;
    }
    if(mech.FBendingType != "" && MParams.FBendingK == 0) {
        cout << "Must set a filament bending constant. Exiting." << endl;
        return false;
    }
    if(mech.FTwistingType != "" && MParams.FTwistingK == 0) {
        cout << "Must set a filament twisting constant. Exiting." << endl;
        return false;
    }
    
    //LINKER
    if(mech.LStretchingType != "" &&
       MParams.LStretchingK.size() != CParams.numLinkerSpecies) {
        cout << "Must set a linker stretching constant for every"
             << " linker species in system. Exiting." << endl;
        return false;
    }
    if(mech.LBendingType != "" &&
       MParams.LBendingK.size() != CParams.numLinkerSpecies) {
        cout << "Must set a linker bending constant for every"
             << " linker species in system. Exiting." << endl;
        return false;
    }
    if(mech.LBendingType != "" &&
       MParams.LBendingTheta.size() != CParams.numLinkerSpecies) {
        cout << "Must set a linker bending angle for every"
             << " linker species in system. Exiting." << endl;
        return false;
    }
    if(mech.LTwistingType != "" &&
       MParams.LTwistingK.size() != CParams.numLinkerSpecies) {
        cout << "Must set a linker twisting constant for every"
        << " linker species in system. Exiting." << endl;
        return false;
    }
    if(mech.LTwistingType != "" &&
       MParams.LTwistingPhi.size() != CParams.numLinkerSpecies) {
        cout << "Must set a linker twisting angle for every"
             << " linker species in system. Exiting." << endl;
        return false;
    }
    
    //MOTOR
    if(mech.MStretchingType != "" &&
       MParams.MStretchingK.size() != CParams.numMotorSpecies) {
        cout << "Must set a motor stretching constant for every"
             << " motor species in system. Exiting." << endl;
        return false;
    }
    if(mech.MBendingType != "" &&
       MParams.MBendingK.size() != CParams.numMotorSpecies) {
        cout << "Must set a motor bending constant for every"
             << " motor species in system. Exiting." << endl;
        return false;
    }
    if(mech.MBendingType != "" &&
       MParams.MBendingTheta.size() != CParams.numMotorSpecies) {
        cout << "Must set a motor bending angle for every"
             << " motor species in system. Exiting." << endl;
        return false;
    }
    if(mech.MTwistingType != "" &&
       MParams.MTwistingK.size() != CParams.numMotorSpecies) {
        cout << "Must set a motor twisting constant for every"
             << " motor species in system. Exiting." << endl;
        return false;
    }
    if(mech.MTwistingType != "" &&
       MParams.MTwistingPhi.size() != CParams.numMotorSpecies) {
        cout << "Must set a motor twisting angle for every"
             << " motor species in system. Exiting." << endl;
        return false;
    }

    
    //BRANCHINGPOINT
    if(mech.BrStretchingType != "" &&
       MParams.BrStretchingK.size() != CParams.numBrancherSpecies) {
        cout << "Must set a brancher stretching constant for every"
             << " brancher species in system. Exiting." << endl;
        return false;
    }
    if(mech.BrBendingType != "" &&
       MParams.BrBendingK.size() != CParams.numBrancherSpecies) {
        cout << "Must set a brancher bending constant for every"
             << " brancher species in system. Exiting." << endl;
        return false;
    }
    if(mech.BrBendingType != "" &&
       MParams.BrBendingTheta.size() != CParams.numBrancherSpecies) {
        cout << "Must set a brancher bending angle for every"
             << " brancher species in system. Exiting." << endl;
        return false;
    }
    if(mech.BrDihedralType != "" &&
       MParams.BrDihedralK.size() != CParams.numBrancherSpecies) {
        cout << "Must set a brancher dihedral constant for every"
             << " brancher species in system. Exiting." << endl;
        return false;
    }
    if(mech.BrPositionType != "" &&
       MParams.BrPositionK.size() != CParams.numBrancherSpecies) {
        cout << "Must set a brancher position constant for every"
             << " brancher species in system. Exiting." << endl;
        return false;
    }
    
    //VOLUME
    if(mech.VolumeFFType != "" && MParams.VolumeK == 0) {
        cout << "Must set a volume constant. Exiting." << endl;
        return false;
    }
    if(mech.VolumeFFType != "" && MParams.VolumeCutoff == 0) {
        cout << "Must set a volume cutoff. Exiting." << endl;
        return false;
    }
    
    //Boundary
    if(mech.BoundaryFFType != "" && BParams.BoundaryK == 0) {
        cout << "Must set a boundary constant. Exiting." << endl;
        return false;
    }
    if(mech.BoundaryFFType != "" && BParams.BScreenLength == 0) {
        cout << "Must set a boundary screen length. Exiting." << endl;
        return false;
    }
    if(mech.BoundaryFFType != "" && BParams.BoundaryCutoff == 0) {
        cout << "Must set a boundary cutoff. Exiting." << endl;
        return false;
    }
    
    return true;
}


bool SystemParameters::checkGeoParameters() {
    
    //Check that grid and compartmentSize match nDim
    if((GParams.nDim == 3 &&
        GParams.NX != 0 && GParams.NY != 0 && GParams.NZ !=0 &&
        GParams.compartmentSizeX != 0 &&
        GParams.compartmentSizeY != 0 &&
        GParams.compartmentSizeZ != 0)){
    }
    else {
        cout << "Grid parameters are invalid. Exiting." << endl;
        return false;
    }
    return true;
}

