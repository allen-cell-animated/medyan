
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_Restart_h
#define MEDYAN_Restart_h

#include "common.h"
#include<stdio.h>
#include "Output.h"
#include "core/controller/MController.h"
#include "core/controller/GController.h"
#include "core/controller/CController.h"
#include "core/controller/DRController.h"
#include <random>
#include <chrono>

#include "core/controller/Controller.h"

#include "Parser.h"
#include "Output.h"
#include "SubSystem.h"
#include "Boundary.h"
#include "CompartmentGrid.h"

#include "FilamentInitializer.h"
#include "BubbleInitializer.h"

#include "Filament.h"
#include "Cylinder.h"
#include "Linker.h"
#include "MotorGhost.h"
#include "BranchingPoint.h"
#include "Bubble.h"
#include "MTOC.h"

#include "SysParams.h"
#include "MathFunctions.h"
#include "core/controller/MController.h"
#include "Cylinder.h"
#include <unordered_set>
#include <unordered_map>
#include 	<tuple>
#include <vector>
#include <algorithm>
using namespace mathfunc;

//FORWARD DECLARATIONS
class SubSystem;
class Cylinder;
class FilamentBindingManager;
class Restart {
private:
    SubSystem *_subSystem; ///< A pointer to the subsystem that this controls
    vector<double> temp_diffrate_vector; ///vector of diffusion rates
    tuple< vector<tuple<short, vector<double>, vector<double>>> , vector<tuple<string, short, vector<vector<double>>>>,
    vector<tuple<string, short, vector<double>>> , vector<vector<double>> > filaments;
    ChemistryData _chemData;
    vector<double> CopyNumbers;
    unordered_multimap<int, tuple<CCylinder*, short>> _unsortedpairings;
    unordered_multimap<int, tuple<CCylinder*, short>> _bunsortedpairings;
    unordered_multimap<int, tuple<int, short>> _singlecylfilunsortedpairings;
    unordered_multimap<int, tuple<int, short>> _bsinglecylfilunsortedpairings;
    vector<LinkerBindingManager*> affectedManagers;
    vector<tuple<string, short, vector<vector<double>>>> boundVector;
    vector<short> branchcylIDs;
    int  _numChemSteps=0;
    
    //gives angle and delta
    vector<double> getAngleDeltaPos(vector<double>leg, vector<double> site1, vector<double> site2){
        vector<double> returnVector;
        double len1, len2, len3;
        vector<double> v1={site1[0]-leg[0],site1[1]-leg[1],site1[2]-leg[2]};
        vector<double> v2={site2[0]-leg[0],site2[1]-leg[1],site2[2]-leg[2]};
        returnVector.push_back(1-acos(std::max(dotProduct(normalizeVector(v1),normalizeVector(v2)),-1.0))*7/22);
        len1=twoPointDistance(leg,site1);
        len2=twoPointDistance(leg,site2);
        len3=twoPointDistance(site1,site2);
        returnVector.push_back(abs(1-(len1+len2)/len3));
        returnVector.push_back(len1);
        returnVector.push_back(len3);
        return returnVector; }
    
    // Goes through single cylinder filaments and decides the appropriate way to activate them.
    void reassignsinglecylfil(bool flag){ //flag 0 - linker/motor, 1-brancher.
        for(auto x:Cylinder::getCylinders()){
            vector<tuple<int, short>> scfmap;
            typedef unordered_multimap<int, tuple<int, short>>:: iterator umit;
            pair<umit, umit> range;
            if(flag ==0)
                range = _singlecylfilunsortedpairings.equal_range(x->getId());
            else
                range = _bsinglecylfilunsortedpairings.equal_range(x->getId());

            vector<int> bVpos; //position in boundVector
            vector<short> bSite; //binding Sites occupied
            vector<short> orginal_bSite; //store orginal bSite
            vector<short> ObSite; //ordered binding Sites occupied
            vector<short> bSitecyl; //all binding sites available in the cylinder.
            vector<short> IDs;
            auto filType = x->getType();
            short deltaBinding = short (SysParams::Geometry().cylinderNumMon[x->getType()] / SysParams::Chemistry().numBindingSites[x->getType()]);
            for(auto it1 = SysParams::Chemistry().bindingSites[filType].begin();
                it1 != SysParams::Chemistry().bindingSites[filType].end(); it1++){
                bSitecyl.push_back((short) *it1);
            }
            for (auto it = range.first; it != range.second; ++it){
                scfmap.push_back(it->second);}//@FOR it
            //TODO if(scfmap.size() > SysParams::Chemistry().numBindingSites[x->getType()])
            //get vectors of bSite
            auto i = 0;
            for(int I=0;I<scfmap.size();I++){
                bVpos.push_back(get<0>(scfmap[I]));
                bSite.push_back(get<1>(scfmap[I]));
                IDs.push_back(i);
                i++;
            }
            //sort in ascending order.
            if(scfmap.size()){
                //if nummonomers is equal to total number of monomers allowed in a cylinder.
                auto nummonomers = min((int) round(x->getMCylinder()->getEqLength()/ SysParams::Geometry().monomerSize[filType]),SysParams::Geometry().cylinderNumMon[filType]);
                if( nummonomers == SysParams::Geometry().cylinderNumMon[filType]){
                    for(auto i = 0; i < bSite.size(); i++){
                        //TODO ensure that proper binding sites are added. Additional binding sites might need to be added as backup.
                        if(flag ==0)
                            _unsortedpairings.insert({bVpos[i],make_tuple(x->getCCylinder(),bSite[i])});
                        else
                            _bunsortedpairings.insert({bVpos[i],make_tuple(x->getCCylinder(),bSite[i])});
                    }
                }
                else{ //in case where there are fewer monomers in the cylinder.
                    auto vecpos = find(branchcylIDs.begin(), branchcylIDs.end(), x->getId());
                    //If the cylinder is both a branch and branching cylinder (like the shaft in letter I), we need to re-do the binding site position when it is branching cylinder.
                    if(vecpos != branchcylIDs.end() && flag == 1){
                          vector<short> posBindingSites=SysParams::Chemistry().bindingSites[filType];
                        //ASSIGN THE BINDING SITE TO BE THE ONE CLOSEST TO WHAT IS GUESSED.
                        for(auto i = 0; i < bSite.size(); i++){
                            int lo=0;
                            int mm;
                            vector<short> test;
                            for(mm=0;mm<posBindingSites.size();mm++){
                                test.push_back(abs(posBindingSites[mm]-bSite[i]));}
                            for(mm=0;mm<posBindingSites.size();mm++){
                                if(test[mm]<test[lo])
                                    lo=mm;}
                            bSite[i] = posBindingSites[lo];
                        }
                    }
                    //Arrange in ascending order.
                    ObSite = bSite;
                    if(bSite.size()>1){
                        short temp;
                        for(i=0;i<ObSite.size();i++){
                            for(auto j=i+1;j<ObSite.size();j++){
                                if(j+1<=ObSite.size()){
                                    if(ObSite[i]>ObSite[j+1])
                                    {
                                        temp = ObSite[i];
                                        ObSite[i] = ObSite[j+1];
                                        ObSite[j+1] = temp;
                                        temp = IDs[j];
                                        IDs[j] = IDs[j+1];
                                        IDs[j+1] = temp;
                                    }
                                }
                            }
                        }
                    }//IF bSite.size() >1
                    //append with other binding sites.
                    while(ObSite.back() < SysParams::Geometry().cylinderNumMon[filType] && ObSite.size() < bSitecyl.size())
                        ObSite.push_back(ObSite.back() + deltaBinding);
                    
                    for(auto i = ObSite.size(); i < bSitecyl.size(); i++){
                        if(ObSite.front() - deltaBinding >0){
                            ObSite.push_back(ObSite.back());
                            for(auto j = ObSite.size()-2; j >0; j--){
                                ObSite.at(j+1)  = ObSite.at(j);
                            }
                            ObSite.at(0)  =ObSite.at(0)- deltaBinding;
                        }
                    }
                    
                    //get COM
                    short mean1 = 0;short mean2 = 0;
                    for(auto i = 0; i < bSitecyl.size(); i++)
                    { mean1 += bSitecyl[i]; mean2 += ObSite[i];}
                    mean1 = mean1/bSitecyl.size();mean2 = mean2/bSitecyl.size();
                    //MAKE SURE THAT YOU ARE NOT MOVING A BRANCH CYLINDER.

                    if(abs(mean1-mean2)!= 0 &&  vecpos != branchcylIDs.end() ){
                        cout<<"Cylinder is not compatible to bind both Link/motor and brancher. Cannot restart. Exiting."<<endl;
                        //exit(EXIT_FAILURE);
                    }
                    else if(abs(mean1-mean2)!=0 && vecpos == branchcylIDs.end()){
                        
                        orginal_bSite = bSite;
                    //move COM to get the necessary translation.
                        for(auto i = 0; i < bSite.size(); i++){
                            bSite[i] = bSite[i] + mean1 - mean2;
                            if(flag ==0){ //if it is a linker/motor
                                //FIND the binding site closest to it.
                                int lo=0;
                                int mm;
                                vector<short> test;
                                for(mm=0;mm<=bSitecyl.size();mm++){
                                    test.push_back(abs(bSitecyl[mm]-bSite[i]));}
                                
                                for(mm=0;mm<=bSitecyl.size();mm++){
                                    if(test[mm]<test[lo])
                                        lo=mm;}

                                _unsortedpairings.insert({bVpos[i],make_tuple(x->getCCylinder(),bSitecyl[lo])});
                            }
                            else
                                //TO DO, check if bSite should be replaced by bSitecyl
                                _bunsortedpairings.insert({bVpos[i],make_tuple(x->getCCylinder(),bSite[i])});
                        }
                        //FIX CCYLINDER
                        auto cc = x->getCCylinder();
                        int nummonomers = min((int) round(x->getMCylinder()->getEqLength()/ SysParams::Geometry().monomerSize[filType]),SysParams::Geometry().cylinderNumMon[filType]);
                        //TURN DOWN OLD MINUS AND PLUS END
                        CMonomer* m1 = cc->getCMonomer(SysParams::Geometry().cylinderNumMon[filType] - nummonomers);
                        m1->speciesMinusEnd(0)->down();
                        m1 = cc->getCMonomer(cc->getSize() - 1);
                        m1->speciesPlusEnd(0)->down();
                        //TURN UP NEW MINUS AND PLUS ENDS.
                        //get the first and last Beads
                        short minus = SysParams::Geometry().cylinderNumMon[filType] - nummonomers + mean1 - mean2;
                        short plus  = SysParams::Geometry().cylinderNumMon[filType] -1 + mean1 - mean2;
                        
                        //check if minus < 0 or plus > cylinder monomer limit
                        //if yes, force minus end and plus end to be in the range
                        if(minus < 0){
                            minus = 0;
                            plus = nummonomers - 1;
                            cout << "Watch out! A minus end index is forced to be 0. Orginal binding site index:" << endl;
                            for (int bindexout = 0; bindexout < orginal_bSite.size(); bindexout++){
                                cout << orginal_bSite[bindexout] << " ";
                            }
                            cout << endl;
                        }
                        else if (plus > (SysParams::Geometry().cylinderNumMon[filType] - 1)){
                            plus = SysParams::Geometry().cylinderNumMon[filType] - 1;
                            minus = SysParams::Geometry().cylinderNumMon[filType] - nummonomers;
                            cout << "Watch out! A plus end index is forced to be " << plus <<". Orginal binding site index:" << endl;
                            for (int bindexout = 0; bindexout < orginal_bSite.size(); bindexout++){
                                cout << orginal_bSite[bindexout] << " ";
                            }
                            cout << endl;
                        }
                        
                        m1 = cc->getCMonomer(minus);
                        m1->speciesMinusEnd(0)->up();
                        m1 = cc->getCMonomer(plus);
                        m1->speciesPlusEnd(0)->up();
                        
                        for(int i = 0; i < cc->getSize(); i++) {
                            if(i>minus && i <plus){ //first CMonomer should be MinusEnd
                                if(cc->getCMonomer(i)->speciesFilament(0)->getN() == 0)
                                    cc->getCMonomer(i)->speciesFilament(0)->up();
                                for(auto j : SysParams::Chemistry().bindingIndices[filType]){
                                    if(cc->getCMonomer(i)->speciesBound(j)->getN() == 0)
                                        cc->getCMonomer(i)->speciesBound(j)->up();}
                            } //@IF
                            else{
                                 if(cc->getCMonomer(i)->speciesFilament(0)->getN() == 1)
                                     cc->getCMonomer(i)->speciesFilament(0)->down();
                                for(auto j : SysParams::Chemistry().bindingIndices[filType]){
                                    if(cc->getCMonomer(i)->speciesBound(j)->getN() == 1)
                                        cc->getCMonomer(i)->speciesBound(j)->down();}
                            } //@ELSE
                        }

                    //FIXED CCYLINDER
                    }
                    else if(abs(mean1-mean2)==0){
                        for(auto i = 0; i < bSite.size(); i++){
                            if(flag ==0)
                                _unsortedpairings.insert({bVpos[i],make_tuple(x->getCCylinder(),bSite[i])});
                            else
                                _bunsortedpairings.insert({bVpos[i],make_tuple(x->getCCylinder(),bSite[i])});
                        }

                    }
                    }// ELSE (IF x->Size > = CYLSIZE)
            }//IF scfmap.size()
            
        } //FOR Cylinders
        
    }//reassign ENDS.
    
    
//cross checks to see linker and motor binding sites have the same distance between them as linker, adds to heap.
    void crosschecklinkermotor(){
        
        short brows=boundVector.size();
        for(int iter=0;iter<=brows-1;iter++){
            vector<tuple<CCylinder*, short>> map;
            auto range = _unsortedpairings.equal_range(iter);
            for (auto it = range.first; it != range.second; ++it){
                map.push_back(it->second);}//@FOR it
            auto b=boundVector.at(iter);
            string boundName=get<0>(b);
            short filamentType=get<1>(b);
            vector<vector<double>> coord=get<2>(b);
            double _numMonPerCyl=SysParams::Geometry().cylinderNumMon[filamentType];
            vector<double> leg1=coord.at(0);
            vector<double> leg2=coord.at(1);
            auto distanceactual=twoPointDistance(leg1,leg2);
            double distanceproj;
            double one,two;
            int check2=0;
            double threshold=0.01;
            for(int I=0;I<map.size();I++){
                for(int J=I+1;J<map.size();J++){
                    auto c1=get<0>(map[I])->getCylinder();
                    auto c2=get<0>(map[J])->getCylinder();
                    auto l1=midPointCoordinate(c1->getFirstBead()->vcoordinate(), c1->getSecondBead()->vcoordinate(),get<1>(map[I])/_numMonPerCyl);
                    auto l2=midPointCoordinate(c2->getFirstBead()->vcoordinate(), c2->getSecondBead()->vcoordinate(),get<1>(map[J])/_numMonPerCyl);
                    distanceproj=twoPointDistance(l1, l2);

                    if(abs((distanceproj-distanceactual)/distanceactual)<threshold)
                    {one=I;two=J;check2=1;threshold=abs((distanceproj-distanceactual)/distanceactual);}
                }}
            if(!check2)
            {   std::cout<<"distances "<<distanceproj<<" "<<distanceactual<<endl;
                cout<<"Serious error! Bound Species (Linker/Motor) with the following coordinates is not bound to a legitimate site. Printing coordinates"<<endl;
                cout<<leg1[0]<<" "<<leg1[1]<<" "<<leg1[2]<<endl;
                cout<<leg2[0]<<" "<<leg2[1]<<" "<<leg2[2]<<endl;
//                exit(EXIT_FAILURE);
            }
            auto c1=get<0>(map[one])->getCylinder();
            auto c2=get<0>(map[two])->getCylinder();
            //@
            if(c1->getId() > c2->getId())        {
                for(auto &Mgr:c1->getCompartment()->getFilamentBindingManagers()){
                    if(dynamic_cast<LinkerBindingManager*>(Mgr.get())) {
                        if(Mgr->getBoundName().compare(boundName)==0){
                            //1. Add necessary diffusing species to the corresponding compartment
                            vector<string> rxnspecies=Mgr->getrxnspecies();
                            setdiffspeciesnumber(rxnspecies,c1);
                            //@
#ifdef NLORIGINAL
                            Mgr->appendpossibleBindings(map[one],map[two]);
#endif
#ifdef NLSTENCILLIST
                            Mgr->appendpossibleBindingsstencil(map[one],map[two]);
#endif
                            _numChemSteps++;}}
                    else if(dynamic_cast<MotorBindingManager*>(Mgr.get())) {
                        if(Mgr->getBoundName().compare(boundName)==0){
                            //2. Add necessary diffusing species to the corresponding compartment
                            vector<string> rxnspecies=Mgr->getrxnspecies();
                            setdiffspeciesnumber(rxnspecies,c1);
                            //@
#ifdef NLORIGINAL
                            Mgr->appendpossibleBindings(map[one],map[two]);
#endif
#ifdef NLSTENCILLIST
                            Mgr->appendpossibleBindingsstencil(map[one],map[two]);
#endif
                            _numChemSteps++;}}
                }//@For Mgr
            }//@IF
            else if(c1->getCompartment() == c2->getCompartment()) {
                for(auto &Mgr:c1->getCompartment()->getFilamentBindingManagers()){
                    if(dynamic_cast<LinkerBindingManager*>(Mgr.get())) {
                        if(Mgr->getBoundName().compare(boundName)==0){
                            //3. Add necessary diffusing species to the corresponding compartment
                            vector<string> rxnspecies=Mgr->getrxnspecies();
                            setdiffspeciesnumber(rxnspecies,c1);
                            //@
#ifdef NLORIGINAL
                            Mgr->appendpossibleBindings(map[two],map[one]);
#endif
#ifdef NLSTENCILLIST
                            Mgr->appendpossibleBindingsstencil(map[two],map[one]);
#endif
                            _numChemSteps++;}}
                    else if(dynamic_cast<MotorBindingManager*>(Mgr.get())) {
                        if(Mgr->getBoundName().compare(boundName)==0){
                            //4. Add necessary diffusing species to the corresponding compartment
                            vector<string> rxnspecies=Mgr->getrxnspecies();
                            setdiffspeciesnumber(rxnspecies,c1);
                            //@
#ifdef NLORIGINAL
                            Mgr->appendpossibleBindings(map[two],map[one]);
#endif
#ifdef NLSTENCILLIST
                            Mgr->appendpossibleBindingsstencil(map[two],map[one]);
#endif
                            _numChemSteps++;}}
                }//@For Mgr
            }//@ELSE IF
            //add in other
            else {
                for(auto &Mgr:c2->getCompartment()->getFilamentBindingManagers()){
                    if(dynamic_cast<LinkerBindingManager*>(Mgr.get())) {
                        if(Mgr->getBoundName().compare(boundName)==0){
                            //5. Add necessary diffusing species to the corresponding compartment
                            vector<string> rxnspecies=Mgr->getrxnspecies();
                            setdiffspeciesnumber(rxnspecies,c2);
                            //@
#ifdef NLORIGINAL
                            Mgr->appendpossibleBindings(map[two],map[one]);
#endif
#ifdef NLSTENCILLIST
                            Mgr->appendpossibleBindingsstencil(map[two],map[one]);
#endif
                            _numChemSteps++;}}
                    else if(dynamic_cast<MotorBindingManager*>(Mgr.get())) {
                        if(Mgr->getBoundName().compare(boundName)==0){
                            //6. Add necessary diffusing species to the corresponding compartment
                            vector<string> rxnspecies=Mgr->getrxnspecies();
                            setdiffspeciesnumber(rxnspecies,c2);
                            //@
#ifdef NLORIGINAL
                            Mgr->appendpossibleBindings(map[two],map[one]);
#endif
#ifdef NLSTENCILLIST
                            Mgr->appendpossibleBindingsstencil(map[two],map[one]);
#endif
                            _numChemSteps++;}}
                }//@For Mgr
            }
    }
    }
    
    //Increases copy number of diffusing species corresponding to bound species in each Compartment by number of events.
    void setdiffspeciesnumber(vector<string> rxnspecies, Cylinder* c){
        int counter=0;
        for(auto sd : _chemData.speciesDiffusing) {
            int events=0;
            for( int it2=0;it2<rxnspecies.size();it2++){
                if((rxnspecies.at(it2)).compare(get<0>(sd))==0){
                    events++;
                }}
            CopyNumbers[counter]=CopyNumbers[counter]-events;
            if(CopyNumbers[counter]<0)
            {cout <<
                "Restart file reaction numbers do not match with diffusing species number."
                << endl;
                exit(EXIT_FAILURE);
                
            }
            events=events+(c->getCompartment()->findSpeciesByName(get<0>(sd)))->getRSpecies().getN();
            (c->getCompartment()->findSpeciesByName(get<0>(sd)))->getRSpecies().setN(events);
            c->getCompartment()->getDiffusionReactionContainer().updatePropensityComprtment();
            counter++;
    }}
public:
    Restart(SubSystem* s, tuple< vector<tuple<short, vector<double>, vector<double>>> , vector<tuple<string, short,
            vector<vector<double>>>> , vector<tuple<string, short, vector<double>>> , vector<vector<double>> > f, ChemistryData _cd)
            : _subSystem(s), filaments(f), _chemData(_cd) {}
    int getnumchemsteps(){return _numChemSteps;}
    void settorestartphase(){
//STEP #1: Get a copy of diffusion rate and set diffusion rate to 0, reset linker motor and branching managers.
    for(auto C : _subSystem->getCompartmentGrid()->getCompartments()) {
        for(auto &it: C->getDiffusionReactionContainer().reactions())
        {temp_diffrate_vector.push_back(it->getRate());
            it->setRate(0.0);}
        C->getDiffusionReactionContainer().updatePropensityComprtment();
        for(auto &Mgr:C->getFilamentBindingManagers()){
#ifdef NLORIGINAL
            Mgr->clearpossibleBindings();
#endif
#ifdef NLSTENCILLIST
            Mgr->clearpossibleBindingsstencil();
#endif
        }}
//STEP #1a: Get cylinders, passivate filament reactions.
//        auto xxx=0;
    for(auto C : _subSystem->getCompartmentGrid()->getCompartments()) {
//        xxx=xxx+C->getCylinders().size();
        for(auto x : C->getCylinders()) {
            x->getCCylinder()->passivatefilreactions();
            x->getCCylinder()->passivatefilcrossreactions();
        }}
//Step #1b. Get copynumber of diffusing species.
    for(auto sd : _chemData.speciesDiffusing) {
        string name = get<0>(sd);
        CopyNumbers.push_back(_subSystem->getCompartmentGrid()->countDiffusingSpecies(name));}
    for(auto C : _subSystem->getCompartmentGrid()->getCompartments()) {
        for(auto sd : _chemData.speciesDiffusing) {
            (C->findSpeciesByName(get<0>(sd)))->getRSpecies().setN(0);}}
//Step #3. Add filament coordinates to be held static during minimization **** NEEEDS TO BE EDITED***
        // coordinates to keep static
        vector<vector<double>> staticbeads=get<3>(filaments);
        int nbeads=get<3>(filaments).size();
        for (int ix=0;ix<nbeads;ix++){
            vector<double> staticbead=staticbeads[ix];
            vector<double> staticbead2;
            staticbead2.insert ( staticbead2.begin(), staticbead.begin()+1, staticbead.begin()+3 );
            for(auto b: Bead::getBeads()) {
                double dis=twoPointDistance(b->vcoordinate(),staticbead2);
                if(dis<=0.00001){
                    b->setstaticstate(true);
                }}}
        
//STEP #2 . updating _possbileBindings of Linkers in each compartment.
        //Filter through probable sites in unsortedpairings by making sure the distance between binding sites
        //in them is the same as bound species bond length
        crosschecklinkermotor();
        crosscheckBranchers();
    }
    
    void addtoHeaplinkermotor(){
//STEP #2. ADD bound Linkers And Motors in inputfile into possible bindings.
        boundVector=get<1>(filaments);
        short brows=boundVector.size();
        vector<vector<double> > site1;
        vector<vector<double> > site2;
        vector<double> angdeltapos;
        for(auto C : _subSystem->getCompartmentGrid()->getCompartments()) {
            for(auto x : C->getCylinders()) {
                vector<double> b1=x->getFirstBead()->vcoordinate();
                vector<double> b2=x->getSecondBead()->vcoordinate();
                // Iterate through each bound coordinate
                for(int iter=0;iter<brows;iter++){
                    auto b=boundVector.at(iter);
                    short filamentType=get<1>(b);
                    if(filamentType==x->getType()){
                        double _numMonPerCyl=0.0;
                        if(x->isMinusEnd() || x->isPlusEnd())
                            _numMonPerCyl=(int) round(x->getMCylinder()->getEqLength()/ SysParams::Geometry().monomerSize[x->getFirstBead()->getType()]);
                        else
                            _numMonPerCyl=SysParams::Geometry().cylinderNumMon[filamentType];
                        vector<short> posBindingSites=SysParams::Chemistry().bindingSites[filamentType];
                        string boundName=get<0>(b);
                        vector<vector<double>> coord=get<2>(b);
                        vector<double> leg1=coord.at(0);
                        vector<double> leg2=coord.at(1);
                        //Leg 1
                        angdeltapos=getAngleDeltaPos(leg1,b1,b2);
                        if( angdeltapos.at(0)<0.001 && angdeltapos.at(1)<0.001){
                            auto f = (Filament*)(x->getParent());

                            double d=0.0;
                            if(_numMonPerCyl< SysParams::Geometry().cylinderNumMon[filamentType]){
                                if(x->isMinusEnd()){
                                    auto vecpos = find(branchcylIDs.begin(), branchcylIDs.end(), x->getId());
                                    if(vecpos!=branchcylIDs.end()) //If it is a branch cylinder, then the CMonomers are re-arranged starting from 0 instead of CMonomer.size().
                                        d = round(angdeltapos.at(2)*SysParams::Geometry().cylinderNumMon[filamentType]/angdeltapos.at(3));
                                        //d = round(angdeltapos.at(2)*_numMonPerCyl/SysParams::Geometry().cylinderSize[filamentType]); //THIS IS THE CORRECT WAY. TEMPORARILY DEPRECATED.
                                    else{
//                                        d = SysParams::Geometry().cylinderNumMon[filamentType] -_numMonPerCyl + round(angdeltapos.at(2)*SysParams::Geometry().cylinderNumMon[filamentType]/angdeltapos.at(3));
                                        d = SysParams::Geometry().cylinderNumMon[filamentType] -_numMonPerCyl + round(angdeltapos.at(2)*_numMonPerCyl/angdeltapos.at(3));
//                                        d = round((1-(angdeltapos.at(3)-angdeltapos.at(2))/SysParams::Geometry().cylinderSize[filamentType])*SysParams::Geometry().cylinderNumMon[filamentType]);
                                    //THIS IS THE CORRECT WAY. TEMPORARILY DEPRECATED.
                                    }
                                    }
                                else
                                    d = round( angdeltapos.at(2)/SysParams::Geometry().monomerSize[x->getFirstBead()->getType()]);
                            }
                            else
                                d = round(angdeltapos.at(2)*_numMonPerCyl/angdeltapos.at(3));
                            if(f->getCylinderVector().size()>1){
                                int lo=0;
                                int mm;
                                vector<short> test;
                                vector<short>bindingsites;
                                bindingsites = posBindingSites;
                                for(mm=0;mm<posBindingSites.size();mm++){
                                    test.push_back(abs(posBindingSites[mm]-d));}
                                //Sort in ascending order.
                                for(mm=0;mm<posBindingSites.size();mm++){
                                    for(int nn=mm+1;nn<posBindingSites.size();nn++){
                                        if(test[nn] < test[mm]){
                                            short temp1, temp2;
                                            temp1 = test[mm];
                                            test[mm] = test[nn];
                                            test[nn] = temp1;
                                            temp2 = bindingsites[mm];
                                            bindingsites[mm] = bindingsites[nn];
                                            bindingsites[nn] = temp2;
                                        }
                                    }
                                }
                                int topn = 2;
                                if(posBindingSites.size()<topn)
                                    topn = posBindingSites.size();
                                for(mm = 0;mm<topn;mm++)
                                    _unsortedpairings.insert({iter,make_tuple(x->getCCylinder(),bindingsites[mm])});
                                //                                for(mm=0;mm<posBindingSites.size();mm++){
                                //                                    if(test[mm]<test[lo])
                                //                                        lo=mm;}
                                //                                _unsortedpairings.insert({iter,make_tuple(x->getCCylinder(),posBindingSites[lo])});
                            }
                            else{
                                _singlecylfilunsortedpairings.insert({x->getId(),make_tuple(iter, d)});
                            }
                        }
                        //@Leg1 ENDS & Leg2
                        angdeltapos=getAngleDeltaPos(leg2,b1,b2);
                        if( angdeltapos.at(0)<0.001 && angdeltapos.at(1)<0.001){
                            auto f = (Filament*)(x->getParent());
                            double d=0.0;
                            if(_numMonPerCyl< SysParams::Geometry().cylinderNumMon[filamentType]){
                                if(x->isMinusEnd()){
                                    auto vecpos = find(branchcylIDs.begin(), branchcylIDs.end(), x->getId());
                                    if(vecpos!=branchcylIDs.end()) //If it is a branch cylinder, then the CMonomers are re-arranged starting from 0 instead of CMonomer.size().
                                        d = round(angdeltapos.at(2)*SysParams::Geometry().cylinderNumMon[filamentType]/angdeltapos.at(3));
                                    //d = round(angdeltapos.at(2)*_numMonPerCyl/SysParams::Geometry().cylinderSize[filamentType]); //THIS IS THE CORRECT WAY. TEMPORARILY DEPRECATED.
                                    else{
//                                        d = SysParams::Geometry().cylinderNumMon[filamentType] -_numMonPerCyl + round(angdeltapos.at(2)*SysParams::Geometry().cylinderNumMon[filamentType]/angdeltapos.at(3));
                          
                                        d = SysParams::Geometry().cylinderNumMon[filamentType] -_numMonPerCyl + round(angdeltapos.at(2)*_numMonPerCyl/angdeltapos.at(3));
//
                                    }
                                    //                                        d = round((1-(angdeltapos.at(3)-angdeltapos.at(2))/SysParams::Geometry().cylinderSize[filamentType])*SysParams::Geometry().cylinderNumMon[filamentType]);
                                    //THIS IS THE CORRECT WAY. TEMPORARILY DEPRECATED.
                                }
                                else
                                    d = round( angdeltapos.at(2)/SysParams::Geometry().monomerSize[x->getFirstBead()->getType()]);
                            }
                            else
                                d = round(angdeltapos.at(2)*_numMonPerCyl/angdeltapos.at(3));
                            if(f->getCylinderVector().size()>1){
                                int lo=0;
                                int mm;
                                vector<short> test;
                                vector<short>bindingsites;
                                bindingsites = posBindingSites;
                                for(mm=0;mm<posBindingSites.size();mm++){
                                    test.push_back(abs(posBindingSites[mm]-d));}
                                //Sort in ascending order.
                                for(mm=0;mm<posBindingSites.size();mm++){
                                    for(int nn=mm+1;nn<posBindingSites.size();nn++){
                                        if(test[nn] < test[mm]){
                                            short temp1, temp2;
                                            temp1 = test[mm];
                                            test[mm] = test[nn];
                                            test[nn] = temp1;
                                            temp2 = bindingsites[mm];
                                            bindingsites[mm] = bindingsites[nn];
                                            bindingsites[nn] = temp2;
                                        }
                                    }
                                }
                                int topn = 2;
                                if(posBindingSites.size()<topn)
                                    topn = posBindingSites.size();
                                for(mm = 0;mm<topn;mm++)
                                    _unsortedpairings.insert({iter,make_tuple(x->getCCylinder(),bindingsites[mm])});
//                                for(mm=0;mm<posBindingSites.size();mm++){
//                                    if(test[mm]<test[lo])
//                                        lo=mm;}
//                                _unsortedpairings.insert({iter,make_tuple(x->getCCylinder(),posBindingSites[lo])});
                            }
                            else{
                                _singlecylfilunsortedpairings.insert({x->getId(),make_tuple(iter, d)});
                            }
                        }//@Leg2 ENDS
                    }//@IF
                }//@for brows
            }}//@Cylinders
        
        //STEP #2 Substep. Check single cylinder filaments to make sure CMonomers are activated appropriately.
        //MinusEnd Cylinders are activated by default with the right most monomer pointing towards the plusEnd.
        reassignsinglecylfil(0);
        }
    
    void addtoHeapbranchers(){
 
        vector<tuple<string, short, vector<double>>> branchVector=get<2>(filaments);
        int iter;
        auto brows=branchVector.size();
        for(auto C : _subSystem->getCompartmentGrid()->getCompartments()) {
            for(auto x : C->getCylinders()) {
                vector<double> b1=x->getFirstBead()->vcoordinate();
                vector<double> b2=x->getSecondBead()->vcoordinate();
                // Iterate through each bound coordinate
                for(iter=0;iter<brows;iter++){
                    auto b=branchVector.at(iter);
                    short filamentType=get<1>(b);
                    if(filamentType==x->getType()){
                        double cylsize=SysParams::Geometry().cylinderSize[filamentType];
                        vector<short> posBindingSites=SysParams::Chemistry().bindingSites[filamentType];
                        double _numMonPerCyl=0.0;
                        if(x->isMinusEnd() || x->isPlusEnd())
                            _numMonPerCyl=(int) (x->getMCylinder()->getEqLength()/SysParams::Geometry().monomerSize[x->getFirstBead()->getType()]);
                        else
                            _numMonPerCyl=SysParams::Geometry().cylinderNumMon[filamentType];
                        auto filamentType = x->getType();
                        string boundName=get<0>(b);
                        vector<double> branch=get<2>(b);
                        //Find the cylinder the brancher is on
                        vector<double> angdeltapos=getAngleDeltaPos(branch,b1,b2);
                        if(angdeltapos.at(0)<0.001 && angdeltapos.at(1)<0.001){
                            auto f = (Filament*)(x->getParent());
                            double d=0.0;
                            if(_numMonPerCyl< SysParams::Geometry().cylinderNumMon[filamentType]){
                                if(x->isMinusEnd())
                                    //TODO MIGHT HAVE TO CHANGE THIS.
                                    d = round((1-(angdeltapos.at(3)-angdeltapos.at(2))/SysParams::Geometry().cylinderSize[filamentType])*SysParams::Geometry().cylinderNumMon[filamentType]);
                                else
                                    d = round( angdeltapos.at(2)/SysParams::Geometry().monomerSize[x->getFirstBead()->getType()]);
                            }
                            else
                                d=round(angdeltapos.at(2)*_numMonPerCyl/angdeltapos.at(3));
                            if(f->getCylinderVector().size()>1){
                                int lo=0;
                                int mm;
                                vector<int> test;
                                for(mm=0;mm<posBindingSites.size();mm++){
                                    test.push_back(abs(posBindingSites[mm]-d));}
                                for(mm=0;mm<posBindingSites.size();mm++){
                                    if(test[mm]<test[lo])
                                        lo=mm;}
                                _bunsortedpairings.insert({iter,make_tuple(x->getCCylinder(),posBindingSites[lo])});
                            }
                            else{

                                _bsinglecylfilunsortedpairings.insert({x->getId(),make_tuple(iter, d)});
                            }
                        }
                        //Find the closest minus end based on distance alone.
                        else if(x->isMinusEnd()&& 0.25>=angdeltapos.at(2)/cylsize){

                            auto f = (Filament*)(x->getParent());
                            
                            if(f->getCylinderVector().size()==1){
                                _bunsortedpairings.insert({iter,make_tuple(x->getCCylinder(),0)});
                                branchcylIDs.push_back(x->getId());
                                
                                auto cc = x->getCCylinder();
                                int nummonomers = min((int) round(x->getMCylinder()->getEqLength()/ SysParams::Geometry().monomerSize[filamentType]),SysParams::Geometry().cylinderNumMon[filamentType]);
                                //TURN DOWN OLD MINUS AND PLUS END
                                CMonomer* m1 = cc->getCMonomer(SysParams::Geometry().cylinderNumMon[filamentType] - nummonomers);
                                m1->speciesMinusEnd(0)->down();
                                m1 = cc->getCMonomer(cc->getSize() - 1);
                                m1->speciesPlusEnd(0)->down();
                                //TURN UP NEW MINUS AND PLUS ENDS.
                                //get the first and last Beads
                                short minus = 0 ;
                                short plus  = nummonomers -1 ;
                                
                                m1 = cc->getCMonomer(minus);
                                if(m1->speciesMinusEnd(0)->getN()!=0)
                                    m1->speciesMinusEnd(0)->down();
                                m1 = cc->getCMonomer(plus);
                                m1->speciesPlusEnd(0)->up();
                            
                                for(int i = 0; i < cc->getSize(); i++) {
                                    if(i>=minus && i <plus){ //first CMonomer should be MinusEnd
                                        if(cc->getCMonomer(i)->speciesFilament(0)->getN() == 0)
                                            cc->getCMonomer(i)->speciesFilament(0)->up();
                                        for(auto j : SysParams::Chemistry().bindingIndices[filamentType]){
                                            if(cc->getCMonomer(i)->speciesBound(j)->getN() == 0)
                                                cc->getCMonomer(i)->speciesBound(j)->up();}
                                    } //@IF
                                    else{
                                        if(cc->getCMonomer(i)->speciesFilament(0)->getN() == 1)
                                            cc->getCMonomer(i)->speciesFilament(0)->down();
                                        for(auto j : SysParams::Chemistry().bindingIndices[filamentType]){
                                            if(cc->getCMonomer(i)->speciesBound(j)->getN() == 1)
                                                cc->getCMonomer(i)->speciesBound(j)->down();}
                                    } //@ELSE
                                }

                            } //IF filament vector has 1 cylinder.
                            else{
                                bool check = false; short sum = 0;
                                int nummonomers = min((int) round(x->getMCylinder()->getEqLength()/ SysParams::Geometry().monomerSize[filamentType]),SysParams::Geometry().cylinderNumMon[filamentType]);
                                for(int i = 0; i < nummonomers; i++) {
                                    sum = sum + x->getCCylinder()->getCMonomer(i)->speciesFilament(0)->getN();
                                }
                                if(x->isMinusEnd() || x->isPlusEnd())
                                    sum++;
                                if(sum == nummonomers )
                                    check = true;
                                if(check){
                                auto m1 = x->getCCylinder()->getCMonomer(0);
                                if(m1->speciesMinusEnd(0)->getN()!=0)
                                    m1->speciesMinusEnd(0)->down();
                                _bunsortedpairings.insert({iter,make_tuple(x->getCCylinder(),0)});
                                }
                                else{
                                    cout<<"A branch filament has more than one cylinder and the minus end is not at CMonomer(0). Cannot restart this file. Exiting."<<endl;
                                    //exit(EXIT_FAILURE);
                                }
                            }
                        }
                    }//@ IF
                }//@ brows
            }//@ Cylinders
        }//@ Compartment
        
        //STEP #2 Substep. Check single cylinder filaments to make sure CMonomers are activated appropriately.
        //MinusEnd Cylinders are activated by default with the right most monomer pointing towards the plusEnd.
        reassignsinglecylfil(1);
    }
    
    void crosscheckBranchers(){
        //Step 3A. Sort through, update possible bindings, fire reaction one by one, handle callback.
        vector<tuple<string, short, vector<double>>> branchVector=get<2>(filaments);
        auto brows=branchVector.size();
        for(auto iter=0;iter<brows;iter++){
            vector<tuple<CCylinder*, short>> map;
            auto range = _bunsortedpairings.equal_range(iter);
            for (auto it = range.first; it != range.second; ++it){
                map.push_back(it->second);}//@FOR it
            auto b=branchVector.at(iter);
            string boundName=get<0>(b);
            short filamentType=get<1>(b);
            vector<double> branch=get<2>(b);
            double _numMonPerCyl=SysParams::Geometry().cylinderNumMon[filamentType];
            double threshold=SysParams::Geometry().cylinderSize[filamentType];
            double one,two;
            int check2=0;
            short pos;
            for(int I=0;I<map.size();I++){
                for(int J=I+1;J<map.size();J++){
                    auto c1=get<0>(map[I])->getCylinder();
                    auto c2=get<0>(map[J])->getCylinder();
                    auto pos1=get<1>(map[I]);
                    auto pos2=get<1>(map[J]);
                    if(c1->getId()!=c2->getId() && c1->getType()==filamentType && c2->getType()==filamentType){
                        vector<double> l1;
                        short check=0;
                        if(pos1!=0)
                        {l1=midPointCoordinate(c1->getFirstBead()->vcoordinate(), c1->getSecondBead()->vcoordinate(),get<1>(map[I])/_numMonPerCyl);check=1;}
                        else if (pos2!=0)
                        {l1=midPointCoordinate(c2->getFirstBead()->vcoordinate(), c2->getSecondBead()->vcoordinate(),get<1>(map[J])/_numMonPerCyl);check=2;}
                        if(check>0){
                            auto distanceproj=twoPointDistance(l1, branch);
                            if(distanceproj>0.0001)
                            {cout<<"Serious error! Brancher "<<iter<<" binding site does not exist"<<endl; break;}
                            if(check==1)
                                distanceproj=twoPointDistance(branch, c2->getFirstBead()->vcoordinate());
                            else if(check==2)
                                distanceproj=twoPointDistance(branch, c1->getFirstBead()->vcoordinate());
                            if(distanceproj<threshold)
                            {if(check==1)
                            {one=I;two=J;pos=pos1;}
                            else
                            {one=J;two=I;pos=pos2;}
                                check2=1;threshold=distanceproj;
                            }}
                    } //@IF get ID
                }}
            if(!check2)
            {cout<<"Serious error! Bound Species (Brancher) with the following coordinates is not bound to a legitimate site"<<endl;
                cout<<branch[0]<<" "<<branch[1]<<" "<<branch[2]<<endl;}
            auto c11=get<0>(map[one])->getCylinder()->getFirstBead()->vcoordinate();
            auto c1=get<0>(map[one])->getCylinder();
            auto c21=get<0>(map[two])->getCylinder()->getFirstBead()->vcoordinate();
            
            for(auto &Mgr:c1->getCompartment()->getFilamentBindingManagers()){
                if(dynamic_cast<BranchingManager*>(Mgr.get())) {
                    if(Mgr->getBoundName().compare(boundName)==0){
                        vector<string> rxnspecies=Mgr->getrxnspecies();
                        int counter=0;
                        for(auto sd : _chemData.speciesDiffusing) {
                            int events=0;
                            for( int it2=0;it2<rxnspecies.size();it2++){
                                if((rxnspecies.at(it2)).compare(get<0>(sd))==0){
                                    events++;
                                }}
                            if((rxnspecies.at(1).compare(get<0>(sd)))!=0)
                                CopyNumbers[counter]=CopyNumbers[counter]-events;
                            events=events+(c1->getCompartment()->findSpeciesByName(get<0>(sd)))->getRSpecies().getN();
                            (c1->getCompartment()->findSpeciesByName(get<0>(sd)))->getRSpecies().setN(events);
                            c1->getCompartment()->getDiffusionReactionContainer().updatePropensityComprtment();
                            counter++;
                        }
#ifdef NLORIGINAL
                        Mgr->appendpossibleBindings(map[one],map[two]);
#endif
#ifdef NLSTENCILLIST
                        Mgr->appendpossibleBindingsstencil(map[one],map[two]);
#endif
                        _numChemSteps++;
                    }}}
        } //@brows
    }
    void redistributediffusingspecies(){
        int numCmprts=(_subSystem->getCompartmentGrid()->getCompartments()).size();
        vector<double> eventVec;
        int counter=0;
        for(int it=0;it<CopyNumbers.size();it++)
            eventVec.push_back(ceil(CopyNumbers[it]/numCmprts));
        for(auto C : _subSystem->getCompartmentGrid()->getCompartments()) {
            counter=0;
            for(auto sd : _chemData.speciesDiffusing) {
                if(CopyNumbers[counter]>0)
                {(C->findSpeciesByName(get<0>(sd)))->getRSpecies().setN(std::min(eventVec[counter],CopyNumbers[counter]));
                    CopyNumbers[counter]=CopyNumbers[counter]-eventVec[counter];}
                else
                    (C->findSpeciesByName(get<0>(sd)))->getRSpecies().setN(0);
                counter++;}}
        counter=0;
        for(auto C : _subSystem->getCompartmentGrid()->getCompartments()){
            for(auto &it: C->getDiffusionReactionContainer().reactions())
            {it->setRate(temp_diffrate_vector[counter]);
                counter++;}
            C->getDiffusionReactionContainer().updatePropensityComprtment();}
    }
};
#endif
