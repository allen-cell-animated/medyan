
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

#ifndef M3SYM_Restart_h
#define M3SYM_Restart_h

#include "common.h"

#include "Output.h"
#include "MController.h"
#include "GController.h"
#include "CController.h"
#include "DRController.h"
#include <random>
#include <chrono>

#include "Controller.h"

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
#include "MController.h"
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
tuple< vector<tuple<short, vector<double>, vector<double>>> , vector<tuple<string, short, vector<vector<double>>>> , vector<tuple<string, short, vector<double>>> , vector<vector<double>> > filaments;
ChemistryData _chemData;
vector<double> CopyNumbers;
unordered_multimap<int, tuple<CCylinder*, short>> _unsortedpairings;
vector<LinkerBindingManager*> affectedManagers;
vector<tuple<string, short, vector<vector<double>>>> boundVector;
int  _numChemSteps=0;
    //gives angle and delta
    vector<double> getAngleDeltaPos(vector<double>leg, vector<double> site1, vector<double> site2){
        vector<double> returnVector;
        double len1, len2, len3;
        vector<double> v1={site1[0]-leg[0],site1[1]-leg[1],site1[2]-leg[2]};
        vector<double> v2={site2[0]-leg[0],site2[1]-leg[1],site2[2]-leg[2]};
        returnVector.push_back(1-acos(std::max(dotProduct(normalizedVector(v1),normalizedVector(v2)),-1.0))*7/22);
        len1=twoPointDistance(leg,site1);
        len2=twoPointDistance(leg,site2);
        len3=twoPointDistance(site1,site2);
        returnVector.push_back(abs(1-(len1+len2)/len3));
        returnVector.push_back(len1);
        returnVector.push_back(len3);
        return returnVector; }
    
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
            double one,two;
            int check2=0;
            double threshold=0.01;
            for(int I=0;I<map.size();I++){
                for(int J=I+1;J<map.size();J++){
                    auto c1=get<0>(map[I])->getCylinder();
                    auto c2=get<0>(map[J])->getCylinder();
                    auto l1=midPointCoordinate(c1->getFirstBead()->coordinate, c1->getSecondBead()->coordinate,get<1>(map[I])/_numMonPerCyl);
                    auto l2=midPointCoordinate(c2->getFirstBead()->coordinate, c2->getSecondBead()->coordinate,get<1>(map[J])/_numMonPerCyl);
                    auto distanceproj=twoPointDistance(l1, l2);
                    if(abs((distanceproj-distanceactual)/distanceactual)<threshold)
                    {one=I;two=J;check2=1;threshold=abs((distanceproj-distanceactual)/distanceactual);}
                }}
            if(!check2)
            {cout<<"Serious error! Bound Species (Linker/Motor) with the following coordinates is not bound to a legitimate site"<<endl;
                cout<<leg1[0]<<" "<<leg1[1]<<" "<<leg1[2]<<endl;
                cout<<leg2[0]<<" "<<leg2[1]<<" "<<leg2[2]<<endl;}
            auto c1=get<0>(map[one])->getCylinder();
            auto c2=get<0>(map[two])->getCylinder();
            //@
            if(c1->getID() > c2->getID())        {
                for(auto &Mgr:c1->getCompartment()->getFilamentBindingManagers()){
                    if(dynamic_cast<LinkerBindingManager*>(Mgr.get())) {
                        if(Mgr->getBoundName().compare(boundName)==0){
                            //1. Add necessary diffusing species to the corresponding compartment
                            vector<string> rxnspecies=Mgr->getrxnspecies();
                            setdiffspeciesnumber(rxnspecies,c1);
                            //@
                            Mgr->appendpossibleBindings(map[one],map[two]);
                            _numChemSteps++;}}
                    else if(dynamic_cast<MotorBindingManager*>(Mgr.get())) {
                        if(Mgr->getBoundName().compare(boundName)==0){
                            //2. Add necessary diffusing species to the corresponding compartment
                            vector<string> rxnspecies=Mgr->getrxnspecies();
                            setdiffspeciesnumber(rxnspecies,c1);
                            //@
                            Mgr->appendpossibleBindings(map[one],map[two]);
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
                            Mgr->appendpossibleBindings(map[two],map[one]);
                            _numChemSteps++;}}
                    else if(dynamic_cast<MotorBindingManager*>(Mgr.get())) {
                        if(Mgr->getBoundName().compare(boundName)==0){
                            //4. Add necessary diffusing species to the corresponding compartment
                            vector<string> rxnspecies=Mgr->getrxnspecies();
                            setdiffspeciesnumber(rxnspecies,c1);
                            //@
                            Mgr->appendpossibleBindings(map[two],map[one]);
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
                            Mgr->appendpossibleBindings(map[two],map[one]);
                            _numChemSteps++;}}
                    else if(dynamic_cast<MotorBindingManager*>(Mgr.get())) {
                        if(Mgr->getBoundName().compare(boundName)==0){
                            //6. Add necessary diffusing species to the corresponding compartment
                            vector<string> rxnspecies=Mgr->getrxnspecies();
                            setdiffspeciesnumber(rxnspecies,c2);
                            //@
                            Mgr->appendpossibleBindings(map[two],map[one]);
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
            events=events+(c->getCompartment()->findSpeciesByName(get<0>(sd)))->getRSpecies().getN();
            (c->getCompartment()->findSpeciesByName(get<0>(sd)))->getRSpecies().setN(events);
//          std::cout<<get<0>(sd)<<" "<<(c1->getCompartment()->findSpeciesByName(get<0>(sd)))->getRSpecies().getN()<<endl;
            c->getCompartment()->getDiffusionReactionContainer().updatePropensityComprtment();
            counter++;
    }}
public:
    Restart(SubSystem* s, tuple< vector<tuple<short, vector<double>, vector<double>>> , vector<tuple<string, short, vector<vector<double>>>> , vector<tuple<string, short, vector<double>>> , vector<vector<double>> > f, ChemistryData _cd) : _subSystem(s), filaments(f), _chemData(_cd) {}
    int getnumchemsteps(){return _numChemSteps;}
    void settorestartphase(){
//STEP #1: Get a copy of diffusion rate and set diffusion rate to 0, reset linker motor and branching managers.
    for(auto C : _subSystem->getCompartmentGrid()->getCompartments()) {
        for(auto &it: C->getDiffusionReactionContainer().reactions())
        {temp_diffrate_vector.push_back(it->getRate());
            it->setRate(0.0);}
        C->getDiffusionReactionContainer().updatePropensityComprtment();
        for(auto &Mgr:C->getFilamentBindingManagers()){Mgr->clearpossibleBindings();
            //            std::cout<<Mgr->numBindingSites()<<endl;
        }}
    
//STEP #1a: Get cylinders, passivate filament reactions.
    for(auto C : _subSystem->getCompartmentGrid()->getCompartments()) {
        for(auto x : C->getCylinders()) {
            x->getCCylinder()->passivatefilreactions();
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
                double dis=twoPointDistance(b->coordinate,staticbead2);
                if(dis<=0.00001){
                    b->setstaticstate(true);
                }}}
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
                vector<double> b1=x->getFirstBead()->coordinate;
                vector<double> b2=x->getSecondBead()->coordinate;
                // Iterate through each bound coordinate
                for(int iter=0;iter<brows;iter++){
                    auto b=boundVector.at(iter);
                    short filamentType=get<1>(b);
                    if(filamentType==x->getType()){
                        vector<short> posBindingSites=SysParams::Chemistry().bindingSites[filamentType];
                        double _numMonPerCyl=SysParams::Geometry().cylinderNumMon[filamentType];
                        string boundName=get<0>(b);
                        vector<vector<double>> coord=get<2>(b);
                        vector<double> leg1=coord.at(0);
                        vector<double> leg2=coord.at(1);
                        //Leg 1
                        angdeltapos=getAngleDeltaPos(leg1,b1,b2);
                        if( angdeltapos.at(0)<0.001 && angdeltapos.at(1)<0.001){
                            double d=round(angdeltapos.at(2)*_numMonPerCyl/angdeltapos.at(3));
                            int lo=0;
                            int mm;
                            vector<short> test;
                            for(mm=0;mm<=posBindingSites.size();mm++){
                                test.push_back(abs(posBindingSites[mm]-d));}
                            for(mm=0;mm<=posBindingSites.size();mm++){
                                if(test[mm]<test[lo])
                                    lo=mm;}
                            _unsortedpairings.insert({iter,make_tuple(x->getCCylinder(),posBindingSites[lo])});}
                        //@Leg1 ENDS & Leg2
                        angdeltapos=getAngleDeltaPos(leg2,b1,b2);
                        if( angdeltapos.at(0)<0.001 && angdeltapos.at(1)<0.001){
                            double d=round(angdeltapos.at(2)*_numMonPerCyl/angdeltapos.at(3));
                            int lo=0;
                            int mm;
                            vector<short> test;
                            for(mm=0;mm<=posBindingSites.size();mm++){
                                test.push_back(abs(posBindingSites[mm]-d));}
                            for(mm=0;mm<=posBindingSites.size();mm++){
                                if(test[mm]<test[lo])
                                    lo=mm;}
                            _unsortedpairings.insert({iter,make_tuple(x->getCCylinder(),posBindingSites[lo])});}
                        //@Leg2 ENDS
                    }//@IF
                }//@for brows
            }}//@Cylinders
        //STEP #2A . updating _possbileBindings of Linkers in each compartment.
        //Filter through probable sites in unsortedpairings by making sure the distance between binding sites in them is the same as bound species bond length
        crosschecklinkermotor();
        }
    
    void addtoHeapbranchers(){
        unordered_multimap<int, tuple<CCylinder*, short>> _bunsortedpairings;
        vector<tuple<string, short, vector<double>>> branchVector=get<2>(filaments);
        int iter;
        auto brows=branchVector.size();
        for(auto C : _subSystem->getCompartmentGrid()->getCompartments()) {
            for(auto x : C->getCylinders()) {
                vector<double> b1=x->getFirstBead()->coordinate;
                vector<double> b2=x->getSecondBead()->coordinate;
                // Iterate through each bound coordinate
                for(iter=0;iter<brows;iter++){
                    auto b=branchVector.at(iter);
                    short filamentType=get<1>(b);
                    if(filamentType==x->getType()){
                        double cylsize=SysParams::Geometry().cylinderSize[filamentType];
                        vector<short> posBindingSites=SysParams::Chemistry().bindingSites[filamentType];
                        double _numMonPerCyl=SysParams::Geometry().cylinderNumMon[filamentType];
                        string boundName=get<0>(b);
                        vector<double> branch=get<2>(b);
                        //Find the cylinder the brancher is on
                        vector<double> angdeltapos=getAngleDeltaPos(branch,b1,b2);
                        if(angdeltapos.at(0)<0.001 && angdeltapos.at(1)<0.001){
                            double d=round(angdeltapos.at(2)*_numMonPerCyl/angdeltapos.at(3));
                            int lo=0;
                            int mm;
                            vector<int> test;
                            for(mm=0;mm<=posBindingSites.size();mm++){
                                test.push_back(abs(posBindingSites[mm]-d));}
                            for(mm=0;mm<=posBindingSites.size();mm++){
                                if(test[mm]==0)
                                    lo=mm;}
                            _bunsortedpairings.insert({iter,make_tuple(x->getCCylinder(),posBindingSites[lo])});}
                        //Find the closest minus end based on distance alone.
                        else if(x->isMinusEnd()&& 0.25>=angdeltapos.at(2)/cylsize){
                            _bunsortedpairings.insert({iter,make_tuple(x->getCCylinder(),0)});
                        }
                    }//@ IF
                }//@ brows
            }//@ Cylinders
        }//@ Compartment
        //Step 3A. Sort through, update possible bindings, fire reaction one by one, handle callback.
        for(iter=0;iter<brows;iter++){
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
                    if(c1->getID()!=c2->getID() && c1->getType()==filamentType && c2->getType()==filamentType){
                        vector<double> l1;
                        short check=0;
                        if(pos1!=0)
                        {l1=midPointCoordinate(c1->getFirstBead()->coordinate, c1->getSecondBead()->coordinate,get<1>(map[I])/_numMonPerCyl);check=1;}
                        else if (pos2!=0)
                        {l1=midPointCoordinate(c2->getFirstBead()->coordinate, c2->getSecondBead()->coordinate,get<1>(map[J])/_numMonPerCyl);check=2;}
                        if(check>0){
                            auto distanceproj=twoPointDistance(l1, branch);
                            if(distanceproj>0.0001)
                            {cout<<"Serious error! Brancher "<<iter<<" binding site does not exist"<<endl; break;}
                            if(check==1)
                                distanceproj=twoPointDistance(branch, c2->getFirstBead()->coordinate);
                            else if(check==2)
                                distanceproj=twoPointDistance(branch, c1->getFirstBead()->coordinate);
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
            auto c11=get<0>(map[one])->getCylinder()->getFirstBead()->coordinate;
            auto c1=get<0>(map[one])->getCylinder();
            auto c21=get<0>(map[two])->getCylinder()->getFirstBead()->coordinate;
            //            std::cout<<c11[0]<<" "<<c11[1]<<" "<<c11[2]<<" & "<<c21[0]<<" "<<c21[1]<<" "<<c21[2]<<" pos "<<p1<<" & "<<p2<<endl;
            for(auto &Mgr:c1->getCompartment()->getFilamentBindingManagers()){
                if(dynamic_cast<BranchingManager*>(Mgr.get())) {
                    if(Mgr->getBoundName().compare(boundName)==0){
                        //std::cout<<Mgr->getBoundName()<<" "<<boundName<<endl;
                        vector<string> rxnspecies=Mgr->getrxnspecies();
                        int counter=0;
                        for(auto sd : _chemData.speciesDiffusing) {
                            int events=0;
                            for( int it2=0;it2<rxnspecies.size();it2++){
                                if((rxnspecies.at(it2)).compare(get<0>(sd))==0){
                                    events++;
                                }}
                            //                            std::cout<<rxnspecies.at(1)<<" "<<get<0>(sd)<<" "<<((rxnspecies.at(1).compare(get<0>(sd))))<<endl;
                            if((rxnspecies.at(1).compare(get<0>(sd)))!=0)
                                CopyNumbers[counter]=CopyNumbers[counter]-events;
                            events=events+(c1->getCompartment()->findSpeciesByName(get<0>(sd)))->getRSpecies().getN();
                            (c1->getCompartment()->findSpeciesByName(get<0>(sd)))->getRSpecies().setN(events);
                            //                        std::cout<<get<0>(sd)<<" "<<(c1->getCompartment()->findSpeciesByName(get<0>(sd)))->getRSpecies().getN()<<endl;
                            c1->getCompartment()->getDiffusionReactionContainer().updatePropensityComprtment();
                            counter++;
                        }
                        Mgr->appendpossibleBindings(map[one],map[two]);
                        //                        _cController->run(1);
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
                //                std::cout<<get<0>(sd)<<" "<<(C->findSpeciesByName(get<0>(sd)))->getRSpecies().getN()<<endl;
                counter++;}}
        counter=0;
        for(auto C : _subSystem->getCompartmentGrid()->getCompartments()){
            for(auto &it: C->getDiffusionReactionContainer().reactions())
            {it->setRate(temp_diffrate_vector[counter]);
                //            std::cout<<it->getRate()<<endl;
                counter++;}
            C->getDiffusionReactionContainer().updatePropensityComprtment();}
    }
};
#endif
