
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



using namespace mathfunc;


//doesn't like to include ReactionContainer.h

class DissipationTracker{
    
private:
    
    SubSystem* _subSystem;   ///< A pointer to subsytem for creation of callbacks, etc.
    ChemistryData _chemData; ///<The chemistry data for the system
    CompartmentGrid* _grid ;
    MController* _mcon;
    
    // reaction counter - not used
    int count;
    
    // cumulative dissipated Gibbs free energy
    double cumDissEnergy;
    
    // Mechanical energy before checmial simulation
    double G1;
    
    // Mechanical energy after chemical simulation, and subsequent mechanical equilibration
    double G2;
    
    // Dissipated chemical free energy
    double GChem;
    

    
    
public:
    
    // Constructor, allow access to objects that information is needed from, set all energy
    // tracking variables to zero
    DissipationTracker(SubSystem* subsystem, ChemistryData chemdata, CompartmentGrid* grid, MController* mcon = NULL): _subSystem(subsystem), _chemData(chemdata), _grid(grid), _mcon(mcon) {
       count = 0;
       cumDissEnergy=0;
       G1=0;
       G2=0;
       GChem=0;
    };

    // increment the reaction counter - not used
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
    
    // Find the Gibbs free energy change for a reaction using its ReactionBase representation
    double getDelGChem(ReactionBase* re){
        
        // get the type of reaction
        ReactionType reType = re->getReactionType();
        
        // get the marker for whether it is rev, irrev, or diff
        RevMType revM = re->getRevMarker();
        
        // get the number of reactants
        int M = re->getM();
        
        // get the number of products
        int N = re->getN();
        
        // calculate sigma
        // int sig = M-N;
        
        // get reaction rate
        float aPlus = re->getBareRate();
        
        // for a vector of stoichiometric coefficients, assumed to be 1 for all
        
        vector<int> reacNu(M,1);
        vector<int> prodNu(N,1);
        
        // get vector of copy numbers of reactants and products
        vector<species_copy_t> reacN = re->getReactantCopyNumbers();
        vector<species_copy_t> prodN = re->getProductCopyNumbers();
        
        vector<string> reacNames = re->getReactantSpecies();
        vector<string> prodNames = re->getProductSpecies();
        
        float delGZero;
        float aMin;
        
        
        // declare delG and set it to 0
        float delG;
        delG=0;
        
        if(reType==0) {
            // Regular Reaction
            if(revM==0){
                // Irreversible
                delGZero =  re->getRevNumber();
                delG = delGIrrChemTherm(delGZero, reacN, reacNu, prodN, prodNu);
                
            } else {
                // Reversible
                aMin = re->getRevNumber();
                delG = delGRevChemTherm(aPlus, aMin, reacN, reacNu,  prodN, prodNu);
        
            }
             
        } else if(reType==1){
             // Diffusion Reaction
             //   delG = delGDifChemTherm(reacN[0],prodN[0]);
             
        } else if(reType==2){
            // Polymerization Plus End
            if(revM==0){
                // Irreversible
                delGZero = re->getRevNumber();
                species_copy_t nMon = reacN[0];
                delG = delGPolyIrrTherm(delGZero,nMon,"P");
                
            } else {
                // Reversible
                aMin = re->getRevNumber();
                species_copy_t nMon;
                nMon = reacN[0];
                delG = delGPolyRevTherm(aPlus,aMin,nMon,"P");
                
            }

        } else if(reType==3){
            // Polymerization Minus End
            if(revM==0){
                // Irreversible
                delGZero = re->getRevNumber();
                species_copy_t nMon = reacN[0];
                delG = delGPolyIrrTherm(delGZero,nMon,"P");
                
            } else {
                // Reversible
                aMin = re->getRevNumber();
                species_copy_t nMon = reacN[0];
                delG = delGPolyRevTherm(aPlus,aMin,nMon,"P");
            }

        } else if(reType==4){
            // Depolymerization Plus End
            if(revM==0){
                // Irreversible
                delGZero = re->getRevNumber();
                species_copy_t nMon = prodN[0];
                delG = delGPolyIrrTherm(delGZero,nMon,"D");
                
            } else {
                // Reversible
                
                aMin = re->getRevNumber();
                species_copy_t nMon = prodN[0];
                delG = delGPolyRevTherm(aPlus,aMin,nMon,"D");
            }
        } else if(reType==5){
            // Depolymerization Minus End
            if(revM==0){
                // Irreversible
                delGZero = re->getRevNumber();
                species_copy_t nMon = prodN[0];
                delG = delGPolyIrrTherm(delGZero,nMon,"D");
                
            } else {
                // Reversible
                aMin = re->getRevNumber();
                species_copy_t nMon = prodN[0];
                delG = delGPolyRevTherm(aPlus,aMin,nMon,"D");
            }
        } else if(reType==6){
            // Linker Binding
            vector<float> vecRates;
            vecRates = re->getLinkerRates();
            aPlus = vecRates[0];
            aMin = vecRates[1];
            
        } else if(reType==7){
            // Motor Binding
            if(revM==0){
                // Irreversible
                
            } else {
                // Reversible
            }
        } else if(reType==8){
            // Linker Unbinding
            
            vector<float> vecRates;
            vecRates = re->getLinkerRates();
            aPlus = vecRates[1];
            aMin = vecRates[0];
            
        } else if(reType==9){
            // Motor Unbinding
            if(revM==0){
                // Irreversible
                
            } else {
                // Reversible
            }
        } else if(reType==10){
            // Motor Walking Forward
            if(revM==0){
                // Irreversible
                
            } else {
                // Reversible
            }
        } else if(reType==11){
            // Motor Walking Backward
            if(revM==0){
                // Irreversible
                
            } else {
                // Reversible
            }
        } else if(reType==12){
            // Filmanet Aging
            if(revM==0){
                // Irreversible
                
                vector<species_copy_t> numR;
                numR.push_back(Filament::countSpecies(0,reacNames[0]));
                
                vector<species_copy_t> numP;
                numP.push_back(Filament::countSpecies(0,prodNames[0]));
                
                delGZero = (re->getRevNumber());
                delG = delGIrrChemTherm(delGZero, numR, reacNu, numP, prodNu);
//                cout<<numR[0]<<endl;
//                cout<<numP[0]<<endl;
                
            } else {
                // Reversible
                aMin = re->getRevNumber();
                delG = delGRevChemTherm(aPlus, aMin, reacN, reacNu,  prodN, prodNu);
                
            }
        } else if(reType==13){
            // Filament Creation
            if(revM==0){
                // Irreversible
                
            } else {
                // Reversible
            }
        } else if(reType==14){
            // Filament Destruction
            if(revM==0){
                // Irreversible
                
            } else {
                // Reversible
            }
        } else if(reType==15){
            // Branching
            if(revM==0){
                // Irreversible
                
            } else {
                // Reversible
            }
        } else if(reType==16){
            // Branch Unbinding
            if(revM==0){
                // Irreversible
                
            } else {
                // Reversible
            }
        } else if(reType==17){
            // Severing
            if(revM==0){
                // Irreversible
                
            } else {
                // Reversible
            }
        }
        
//        if(delG>0.01){
//            if(reType!=1){
//                cout<<reType<<endl;
//                cout<<delG<<endl;
//                cout<<GChem<<endl;
//            }
//        }
        
        return delG;
        
        
    }
    
        
    
    
    // increment the GChem counter when a reaction fires
    void updateDelGChem(ReactionBase* re){
        GChem += getDelGChem(re);
    }

    // return the mechanical energy of the system
    double getMechEnergy(){
        return _mcon->getEnergy()/kT;
    }
    
    // return the cumulative dissipated Gibbs free energy
    double getEnergy(){
        return cumDissEnergy;
    }
    
    double getGChemEn(){
        return GChem;
    }
    
    // set the value of G1
    void setG1(){
        G1=getMechEnergy()/kT;

    }
    
    // set the value of G2
    void setG2(){
        G2=getMechEnergy()/kT;
        
    }
    
    // increment cumDissEn after an iteration step has occured
    void updateCumDissEn(){
        cumDissEnergy += GChem + G2 - G1;
    }
    
    // set new values of energy trackers after an iteration step has occured
    void resetAfterStep(){
        GChem=0;
        G1=G2;
        G2=0;
    };

    
    
    
};


#endif
