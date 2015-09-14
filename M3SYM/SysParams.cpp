
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

#include "SysParams.h"

bool SysParams::checkChemParameters(ChemistryData& chem) {
    
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
        return false;
    }
    if(CParams.numFilamentSpecies == 0) {
        cout << "Must provide at least one filament species. Exiting."
        << endl;
        return false;
        
    }
    if(CParams.numPlusEndSpecies == 0) {
        cout << "Must provide at least one plus end species. Exiting."
        << endl;
        return false;
        
    }
    if(CParams.numMinusEndSpecies == 0) {
        cout << "Must provide at least one minus end species. Exiting."
        << endl;
        return false;
    }
    //at least one bound
    if(CParams.numBoundSpecies < 1) {
        cout << "Must provide at least one bound species, which is an empty bounding site. Exiting." << endl;
        return false;
    }
    
    //additional motor params
    if(chem.speciesMotor.size() != CParams.motorNumHeadsMin.size()) {
        
        cout << "Number of minimum motor heads in chemistry input does not "
             << "match the number of motor species. Check these parameters. Exiting."
        <<endl;
        
        return false;
    }
    if(chem.speciesMotor.size() != CParams.motorNumHeadsMax.size()) {
        
        cout << "Number of maximum motor heads in chemistry input does not "
             << "match the number of motor species. Check these parameters. Exiting."
        <<endl;
        return false;
    }
    if(chem.speciesMotor.size() != CParams.motorStepSize.size()) {
        
        cout << "Number of motor step sizes in chemistry input does not "
             << "match the number of motor species. Check these parameters. Exiting."
        <<endl;
        return false;
    }
    
    //check if binding sites are valid
    if(chem.speciesBrancher.size() > 0 && chem.brancherBindingSite == "") {
        cout << "A brancher binding site must be set for brancher species. Exiting."
             << endl;
        return false;
    }
    
    if(chem.speciesLinker.size() > 0 && chem.linkerBindingSite == "") {
        cout << "A linker binding site must be set for linker species. Exiting."
        << endl;
        return false;
    }
    
    if(chem.speciesMotor.size() > 0 && chem.motorBindingSite == "") {
        cout << "A motor binding site must be set for motor species. Exiting."
        << endl;
        return false;
    }
    
    return true;
}

bool SysParams::checkMechParameters(MechanicsFFType& mech) {
    
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
        cout << "Number of linker stretching constants does not match the number of"
             << " linker species in system. Exiting." << endl;
        return false;
    }
    if(mech.LBendingType != "" &&
       MParams.LBendingK.size() != CParams.numLinkerSpecies) {
        cout << "Number of linker bending constants does not match the number of"
             << " linker species in system. Exiting." << endl;
        return false;
    }
    if(mech.LBendingType != "" &&
       MParams.LBendingTheta.size() != CParams.numLinkerSpecies) {
        cout << "Number of linker bending angles does not match the number of"
             << " linker species in system. Exiting." << endl;
        return false;
    }
    if(mech.LTwistingType != "" &&
       MParams.LTwistingK.size() != CParams.numLinkerSpecies) {
        cout << "Number of linker twisting constants does not match the number of"
        << " linker species in system. Exiting." << endl;
        return false;
    }
    if(mech.LTwistingType != "" &&
       MParams.LTwistingPhi.size() != CParams.numLinkerSpecies) {
        cout << "Number of linker twisting angles does not match the number of"
             << " linker species in system. Exiting." << endl;
        return false;
    }
    
    //MOTOR
    if(mech.MStretchingType != "" &&
       MParams.MStretchingK.size() != CParams.numMotorSpecies) {
        cout << "Number of motor stretching constants does not match the number of"
             << " motor species in system. Exiting." << endl;
        return false;
    }
    if(mech.MBendingType != "" &&
       MParams.MBendingK.size() != CParams.numMotorSpecies) {
        cout << "Number of motor bending constants does not match the number of"
             << " motor species in system. Exiting." << endl;
        return false;
    }
    if(mech.MBendingType != "" &&
       MParams.MBendingTheta.size() != CParams.numMotorSpecies) {
        cout << "Number of motor bending angles does not match the number of"
             << " motor species in system. Exiting." << endl;
        return false;
    }
    if(mech.MTwistingType != "" &&
       MParams.MTwistingK.size() != CParams.numMotorSpecies) {
        cout << "Number of motor twisting constants does not match the number of"
             << " motor species in system. Exiting." << endl;
        return false;
    }
    if(mech.MTwistingType != "" &&
       MParams.MTwistingPhi.size() != CParams.numMotorSpecies) {
        cout << "Number of motor twisting angles does not match the number of"
             << " motor species in system. Exiting." << endl;
        return false;
    }

    
    //BRANCHINGPOINT
    if(mech.BrStretchingType != "" &&
       MParams.BrStretchingK.size() != CParams.numBrancherSpecies) {
        cout << "Number of branching point stretching constants does not match the number of"
             << " brancher species in system. Exiting." << endl;
        return false;
    }
    if(mech.BrStretchingType != "" &&
       MParams.BrStretchingL.size() != CParams.numBrancherSpecies) {
        cout << "Number of branching point stretching length does not match the number of"
             << " brancher species in system. Exiting." << endl;
        return false;
    }
    if(mech.BrBendingType != "" &&
       MParams.BrBendingK.size() != CParams.numBrancherSpecies) {
        cout << "Number of branching point bending constants does not match the number of"
             << " brancher species in system. Exiting." << endl;
        return false;
    }
    if(mech.BrBendingType != "" &&
       MParams.BrBendingTheta.size() != CParams.numBrancherSpecies) {
        cout << "Number of branching point bending angles does not match the number of"
             << " brancher species in system. Exiting." << endl;
        return false;
    }
    if(mech.BrDihedralType != "" &&
       MParams.BrDihedralK.size() != CParams.numBrancherSpecies) {
        cout << "Number of branching point dihedral constants does not match the number of"
             << " brancher species in system. Exiting." << endl;
        return false;
    }
    if(mech.BrPositionType != "" &&
       MParams.BrPositionK.size() != CParams.numBrancherSpecies) {
        cout << "Number of branching point position constants does not match the number of"
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

bool SysParams::checkGeoParameters() {
    
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

bool SysParams::checkDyRateParameters(DynamicRateTypes& dy) {
    
    //check types match number of species
    if(dy.dLUnbindingType.size() != CParams.numLinkerSpecies &&
       !dy.dLUnbindingType.empty()) {
        cout << "Number of linker dynamic rate unbinding forms must" <<
                " match the number of species. Exiting." << endl;
        return false;
    }
    
    if(dy.dMUnbindingType.size() != CParams.numMotorSpecies &&
       !dy.dMUnbindingType.empty()) {
        cout << "Number of motor dynamic rate unbinding forms must" <<
                " match the number of species. Exiting." << endl;
        return false;
    }
    if(dy.dMWalkingType.size() != CParams.numMotorSpecies &&
       !dy.dMWalkingType.empty()) {
        cout << "Number of motor dynamic rate walking forms must" <<
                " match the number of species. Exiting." << endl;
        return false;
    }
    
    //now check parameters
    if(dy.dFPolymerizationType != "" && DRParams.dFilPolymerizationCharLength == 0) {
        cout << "Must set a dynamic rate polymerization length. Exiting." << endl;
        return false;
    }
    
    auto numCharLengths = 0;
    auto numAmps = 0;
    
    for(auto &changer : dy.dLUnbindingType) {
        
        if(changer == "BASICCATCHSLIP") {
            numCharLengths += 2;
            numAmps += 2;
        }
        else if(changer == "BASICSLIP") {
            numCharLengths += 1;
        }
        
    }
    if(numCharLengths != SysParams::DynamicRates().
                         dLinkerUnbindingCharLength.size()) {
        cout << "Number of characteristic lengths specified for chosen "
             << "linker unbinding dynamic rate forms is not accurate. Exiting."
        << endl;
        return false;
    }
    
    if(numAmps != SysParams::DynamicRates().
                  dLinkerUnbindingAmplitude.size()) {
        
        
        cout << "Number of amplitudes specified for chosen "
             << "linker unbinding dynamic rate forms is not accurate. Exiting."
        << endl;
        return false;
    }
    
    if(dy.dMUnbindingType.size() != SysParams::DynamicRates().
                                    dMotorUnbindingCharForce.size()) {
        cout << "Number of characteristic forces specified for chosen "
             << "motor unbinding dynamic rate forms is not accurate. Exiting."
        << endl;
        return false;
    }
    
    if(dy.dMWalkingType.size() != SysParams::DynamicRates().
                                  dMotorWalkingCharForce.size()) {
        cout << "Number of characteristic forces specified for chosen "
             << "motor walking dynamic rate forms is not accurate. Exiting."
        << endl;
        return false;
    }
    
    return true;
}

MechParams   SysParams::MParams;
ChemParams   SysParams::CParams;
GeoParams    SysParams::GParams;
BoundParams  SysParams::BParams;
DyRateParams SysParams::DRParams;


