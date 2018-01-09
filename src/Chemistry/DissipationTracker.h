
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_DissipationTracker_h
#define MEDYAN_DissipationTracker_h

#include "common.h"

#include "MathFunctions.h"
#include "ChemSim.h"
#include "ReactionBase.h"
#include "CompartmentGrid.h"
#include "Filament.h"
#include "MController.h"
#include "common.h"
#include "SysParams.h"
#include <fstream>
#include "MotorGhost.h"



using namespace mathfunc;


//doesn't like to include ReactionContainer.h

class DissipationTracker{
    
private:
    

    MController* _mcon;
    
    // reaction counter - not used
    int count;
    
    // cumulative dissipated Gibbs free energy
    double cumDissEnergy;
    
    // Mechanical energy before checmial simulation
    double G1;
    
    // Mechanical energy after chemical simulation
    double GMid;
    
    // Mechanical energy after subsequent mechanical equilibration
    double G2;
    
    // Change in chemical free energy after chemical simulation
    double GChem;
    
    // cumulative dissiapted chemical energy
    double cumDissChemEnergy;
    
    // cumulative dissiapted mechanical energy
    double cumDissMechEnergy;

    // cumulative change in chemical energy
    double cumGChemEn;
    
    // cumulative change in mechanical energy
    double cumGMechEn;
    
    
public:
    
    // Constructor, allow access to objects that information is needed from, set all energy
    // tracking variables to zero
    DissipationTracker(MController* mcon = NULL):  _mcon(mcon) {
       count = 0;
       cumDissEnergy=0;
       cumDissChemEnergy=0;
       cumDissMechEnergy=0;
       cumGChemEn=0;
       cumGMechEn=0;
       G1=0;
       G2=0;
       GMid=0;
       GChem=0;
    };

    // increment the reaction counter
    void updateCount(){
        count++;
    };
    
    // return reaction counter
    int getCount(){return count;};
    
    // print reaction counter
    void printCount(ofstream& outputFile){
        outputFile << count << endl;};
    
    // get the number of compartments in each dimension to form the ratio v used later
    int nx = SysParams::Geometry().NX;
    int ny = SysParams::Geometry().NY;
    int nz = SysParams::Geometry().NZ;
    float v = (float)(nx*ny*nz);
    
    // get the stepFrac parameter for myosin walking, assume all fil and motor are type 0
    double d_step = SysParams::Chemistry().motorStepSize[0];
    double d_total = (double)SysParams::Geometry().cylinderSize[0] /
    SysParams::Chemistry().numBindingSites[0];
    double _stepFrac = d_step / d_total;
    
    // Find the Gibbs free energy change for a reaction using its ReactionBase representation
    double getDelGChem(ReactionBase* re){
        
        // get the type of reaction
        ReactionType reType = re->getReactionType();
        
        
        // get the number of reactants
        int M = re->getM();
        
        // get the number of products
        int N = re->getN();
 
        
        
        // for a vector of stoichiometric coefficients, assumed to be 1 for all
        
        vector<int> reacNu(M,1);
        vector<int> prodNu(N,1);
        
        // get vector of copy numbers of reactants and products
        vector<species_copy_t> reacN = re->getReactantCopyNumbers();
        vector<species_copy_t> prodN = re->getProductCopyNumbers();
        
        vector<string> reacNames = re->getReactantSpecies();
        vector<string> prodNames = re->getProductSpecies();
        
        float delGZero;

        
        // declare delG and set it to 0
        float delG;
        delG=0;
        
        if(reType==0) {
            // Regular Reaction
            delGZero =  re->getGNumber();
            delG = delGGenChem(delGZero, reacN, reacNu, prodN, prodNu);
             
        } else if(reType==1){
             // Diffusion Reaction
                delG = delGDifChem(reacN[0],prodN[0]);
             
        } else if(reType==2){
            // Polymerization Plus End

            delGZero = re->getGNumber();
            species_copy_t nMon = reacN[0];
            delG = delGPolyChem(delGZero,nMon,"P");


        } else if(reType==3){
            // Polymerization Minus End

            delGZero = re->getGNumber();
            species_copy_t nMon = reacN[0];
            delG = delGPolyChem(delGZero,nMon,"P");


        } else if(reType==4){
            // Depolymerization Plus End

            delGZero = re->getGNumber();
            species_copy_t nMon = prodN[0];
            delG = delGPolyChem(delGZero,nMon,"D");

        } else if(reType==5){
            // Depolymerization Minus End

            delGZero = re->getGNumber();
            species_copy_t nMon = prodN[0];
            delG = delGPolyChem(delGZero,nMon,"D");

        } else if(reType==6){
            // Linker Binding
            delGZero=re->getGNumber();
            species_copy_t nMon = reacN[1];
            delG = delGPolyChem(delGZero,nMon,"P");
            
            
        } else if(reType==7){
            // Motor Binding
            double rn=re->getGNumber();
 
            double nh1 = SysParams::Chemistry().motorNumHeadsMin[0];
            double nh2 = SysParams::Chemistry().motorNumHeadsMax[0];
            double nh = (nh1+nh2)/2.0;
            
            delG = delGMyoChem(nh,rn);
            
            
        } else if(reType==8){
            // Linker Unbinding
            delGZero=re->getGNumber();
            species_copy_t nMon = prodN[0];
            delG = delGPolyChem(delGZero,nMon,"D");
            

            
        } else if(reType==9){
            // Motor Unbinding
            double rn=re->getGNumber();
            
            CBound* CBound = re->getCBound();
            SpeciesBound* sm1 = CBound->getFirstSpecies();
            MotorGhost* m = ((CMotorGhost*)sm1->getCBound())->getMotorGhost();
            int nh = m->getNumHeads();
            
            delG = delGMyoChem(nh,rn);
            delG = -delG;

            
        } else if(reType==10){
            // Motor Walking Forward
            delG = re->getGNumber();
            delG = delG*(1/_stepFrac);
            
        } else if(reType==11){
            // Motor Walking Backward
            delG = -(re->getGNumber());
            delG = delG*(1/_stepFrac);
            
        } else if(reType==12){
            // Filament Aging
 
                vector<species_copy_t> numR;
                numR.push_back(Filament::countSpecies(0,reacNames[0]));
            
                vector<species_copy_t> numP;
                numP.push_back(Filament::countSpecies(0,prodNames[0]));
                
                delGZero = (re->getGNumber());
                delG = delGGenChem(delGZero, numR, reacNu, numP, prodNu);
            
        } else if(reType==13){
            // Filament Creation
            
        } else if(reType==14){
            // Filament Destruction
            
        } else if(reType==15){
            // Branching
            
        } else if(reType==16){
            // Branch Unbinding
            
        } else if(reType==17){
            // Severing
            
        }
        
        // check for nan
        
        if(delG<-pow(10,7)){
            cout<<"neg inf"<<endl;
            cout<<"retype is "<<reType<<endl;
            delG=0;
        }
        
        if(delG>pow(10,7)){
            cout<<"pos inf"<<endl;
            cout<<"retype is "<<reType<<endl;
            delG=0;
        }
        
        if(delG!=delG){
            cout<<"nan happened"<<endl;
            delG=0;
        }
        
        return delG;
       
        
    }
    
        
    
    
    // increment the GChem counter when a reaction fires
    void updateDelGChem(ReactionBase* re){
        GChem += getDelGChem(re);
        updateCount();

    }

    // return the mechanical energy of the system
    double getMechEnergy(){
        double ret = _mcon->getEnergy();
        return float(ret/kT);
    }
    
    // return the cumulative dissipated Gibbs free energy
    double getCumDissEnergy(){
        return cumDissEnergy;
    }
    
    // return the cumulative dissipated chemical Gibbs free energy
    double getCumDissChemEnergy(){
        return cumDissChemEnergy;
    }
    
    // return the cumulative dissipated mechanical Gibbs free energy
    double getCumDissMechEnergy(){
        return cumDissMechEnergy;
    }
    
    double getCumGChemEn(){
        return cumGChemEn;
    }
    
    double getCumGMechEn(){
        return cumGMechEn;
    }
    
    // set the value of G1
    void setG1(){
        G1=getMechEnergy();

    }
    
    // set the value of G2
    void setG2(){
        G2=getMechEnergy();
//        cout<<"G1 is "<<G1<<endl;
//        cout<<"G2 is "<<G2<<endl;
        
    }
    
    void setGMid(){
        GMid=getMechEnergy();
    }
    
    // increment cumDissEn after an iteration step has occured
    void updateCumDissEn(){
        cumDissEnergy += GChem + G2 - G1;
    }
    
    // increment cumDissChemEnergy
    void updateCumDissChemEnergy(){
        cumDissChemEnergy += GChem - (G1-GMid);
    }
    
    // increment cumDissMechEnergy
    void updateCumDissMechEnergy(){
        cumDissMechEnergy += G2-GMid;
    }
    
    // increment cumDissMechEnergy
    void updateCumGChemEn(){
        cumGChemEn += GChem;
    }
    
    // increment cumDissMechEnergy
    void updateCumGMechEn(){
        cumGMechEn += G2-G1;
    }
    
    // set new values of energy trackers after an iteration step has occured
    void resetAfterStep(){
        GChem=0;
        GMid=0;
        G1=G2;
    };

    
    
    
};


#endif
