
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
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
#include "MController.h"
#include "GController.h"
#include "CController.h"
#include "DRController.h"
#include <random>
#include <chrono>

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
#include "RestartParams.h"
#include <unordered_set>
#include <unordered_map>
#include 	<tuple>
#include <algorithm>
using namespace mathfunc;

//FORWARD DECLARATIONS
class SubSystem;
class Cylinder;
class FilamentBindingManager;

class Restart {
private:
    SubSystem *_subSystem; ///< A pointer to the subsystem that this controls
    vector<floatingpoint> temp_diffrate_vector; ///vector of diffusion rates
    tuple< vector<tuple<short, vector<floatingpoint>, vector<floatingpoint>>> , vector<tuple<string, short, vector<vector<floatingpoint>>>>,
    vector<tuple<string, short, vector<floatingpoint>>> , vector<vector<floatingpoint>> > filaments;
    ChemistryData _chemData;
    vector<floatingpoint> CopyNumbers;
    unordered_multimap<int, tuple<CCylinder*, short>> _unsortedpairings;
    unordered_multimap<int, tuple<CCylinder*, short>> _bunsortedpairings;
    unordered_multimap<int, tuple<int, short>> _singlecylfilunsortedpairings;
    unordered_multimap<int, tuple<int, short>> _bsinglecylfilunsortedpairings;
    vector<LinkerBindingManager*> affectedManagers;
    vector<tuple<string, short, vector<vector<floatingpoint>>>> boundVector;
    vector<short> branchcylIDs;
    int  _numChemSteps=0;

	fstream _inputFile; ///< input file being used
	restartSystemData _rsystemdata;
	restartBeadData _rBData;
	vector<restartFilData> _rFDatavec;
    vector<restartCylData> _rCDatavec;
    vector<restartMotorData> _rMDatavec;
    vector<restartLinkerData> _rLDatavec;
    vector<restartBrancherData> _rBDatavec;
    vector<restartCompartmentDiffusingData> _rCmpDDatavec;
	restartBulkData _rbdata;
    restartcopynumberrallyData _rtallydata;
    //gives angle and delta
    vector<floatingpoint> getAngleDeltaPos(vector<floatingpoint>leg, vector<floatingpoint> site1, vector<floatingpoint> site2){
        vector<floatingpoint> returnVector;
        floatingpoint len1, len2, len3;
        vector<floatingpoint> v1={site1[0]-leg[0],site1[1]-leg[1],site1[2]-leg[2]};
        vector<floatingpoint> v2={site2[0]-leg[0],site2[1]-leg[1],site2[2]-leg[2]};
        returnVector.push_back(1-acos(std::max<floatingpoint>(dotProduct(normalizeVector
        (v1), normalizeVector(v2)),(floatingpoint)-1.0))*7/22);
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

    //Increases copy number of diffusing species corresponding to bound species in each Compartment by number of events.
    void setdiffspeciesnumber(string diffusingspeciesname, Cylinder* c){
        int counter=0;
        for(auto sd : _chemData.speciesDiffusing) {
            int events=0;
            if(diffusingspeciesname.compare(get<0>(sd))==0){
            	events++;
            }
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
    }
    }
public:
    Restart(SubSystem* s, ChemistryData _cd, const string inputFileName)
            : _subSystem(s), _chemData(_cd) {
	    _inputFile.open(inputFileName);
	    if(!_inputFile.is_open()) {
		    cout << "There was an error parsing file " << inputFileName
		         << ". Exiting." << endl;
		    exit(EXIT_FAILURE);
	    }

    }

    ~Restart(){_inputFile.close();}

    void readNetworkSetup();

    void setupInitialNetwork();

    floatingpoint getrestartime(){return _rsystemdata.time;}

    void settorestartphase(){
//STEP #1: Get a copy of diffusion rate and set diffusion rate to 0, reset linker motor and branching managers.
        for(auto C : _subSystem->getCompartmentGrid()->getCompartments()) {
            for(auto &it: C->getDiffusionReactionContainer().reactions())
            {temp_diffrate_vector.push_back(it->getRate());
                it->setRateMulFactor(0.0f, ReactionBase::RESTARTPHASESWITCH);}
            C->getDiffusionReactionContainer().updatePropensityComprtment();
            for(auto &Mgr:C->getFilamentBindingManagers()){
#ifdef NLORIGINAL
                Mgr->clearpossibleBindings();
#else
                Mgr->clearpossibleBindingsstencil();
#endif
            }}
//STEP #1a: Get cylinders, passivate filament reactions.
        for(auto C : _subSystem->getCompartmentGrid()->getCompartments()) {
            for(auto x : C->getCylinders()) {
                x->getCCylinder()->passivatefilreactions();
                x->getCCylinder()->passivatefilcrossreactions();
            }}
//Step #1b. Get copynumber of diffusing species. This is used later for book keeping
// purposes.
        for(auto sd : _chemData.speciesDiffusing) {
            string name = get<0>(sd);
            CopyNumbers.push_back(_subSystem->getCompartmentGrid()->countDiffusingSpecies(name));}
        //Set copy number of diffusing species in each compartment to 0.
        for(auto C : _subSystem->getCompartmentGrid()->getCompartments()) {
            for(auto sd : _chemData.speciesDiffusing) {
                (C->findSpeciesByName(get<0>(sd)))->getRSpecies().setN(0);}
            for(auto sd : _chemData.speciesBulk) {
                (C->findSpeciesByName(get<0>(sd)))->getRSpecies().setN(0);}
        }
//Step #3. Add filament coordinates to be held static during minimization **** NEEEDS TO BE EDITED***
        // coordinates to keep static
        vector<vector<floatingpoint>> staticbeads=get<3>(filaments);
        int nbeads=get<3>(filaments).size();
        for (int ix=0;ix<nbeads;ix++){
            vector<floatingpoint> staticbead=staticbeads[ix];
            vector<floatingpoint> staticbead2;
            staticbead2.insert ( staticbead2.begin(), staticbead.begin()+1, staticbead.begin()+3 );
            for(auto b: Bead::getBeads()) {
                floatingpoint dis=twoPointDistance(b->vcoordinate(),staticbead2);
                if(dis<=0.00001){
                    b->setstaticstate(true);
                }}}

//STEP #2 . updating _possbileBindings of Linkers in each compartment.
        //Filter through probable sites in unsortedpairings by making sure the distance between binding sites
        //in them is the same as bound species bond length
        addtoHeapLinkerMotorBrancher();
    }

    int getnumchemsteps(){return _numChemSteps;}

    void addtoHeapLinkerMotorBrancher();

    void CBoundinitializerestart();

    void restartupdateCopyNumbers(){
    	// If user requests to use restart file copy numbers after restart.
    	if(SysParams::RUNSTATE == false && SysParams::USECHEMCOPYNUM == false) {
		    //Loop through restart data
		    // Find corresponding compartment and add species to it.
		    for (auto &Cmp:_subSystem->getCompartmentGrid()->getCompartments()) {
			    auto cdiffdata = _rCmpDDatavec[Cmp->getId()];
			    auto diffspeciesnamevec = cdiffdata.speciesnamevec;
			    auto diffspeciescpyvec = cdiffdata.copynumvec;
			    if (Cmp->getId() == cdiffdata.id) {
				    for (auto count = 0; count < diffspeciesnamevec.size(); count++) {
					    Species *sp = Cmp->findSpeciesByName(diffspeciesnamevec[count]);
					    if (sp != nullptr) {
						    sp->getRSpecies().setN(diffspeciescpyvec[count]);
					    } else {
						    LOG(ERROR) << "Cannot find diffusing species of name "
						               << diffspeciesnamevec[count]
						               << " in Compartment ID " << Cmp->getId()
						               << ". Check chemistry input file. Exiting" << endl;
						    throw std::logic_error("Restart unsuccessful!");
					    }
				    }
			    }
		    }
		    int counter = 0;
		    for (auto C : _subSystem->getCompartmentGrid()->getCompartments()) {
			    for (auto &it: C->getDiffusionReactionContainer().reactions()) {
				    it->setRateMulFactor(1.0f, ReactionBase::RESTARTPHASESWITCH);
			    }
			    C->getDiffusionReactionContainer().updatePropensityComprtment();
		    }
		    //Bulk Species
		    for (auto &s : _chemData.speciesBulk) {

		    }
	    }
	    checktallyofcopynumbers();
    }

    void checktallyofcopynumbers(){
        for(int i =0; i < _rtallydata.copynumvec.size(); i++){
            string name = _rtallydata.speciesnamevec.at(i);
            auto copyNum = _subSystem->getCompartmentGrid()->countDiffusingSpecies(name);
            if(_rtallydata.copynumvec[i] != copyNum){
                LOG(ERROR) << "Species copy number does not tally. Species name "
                           << name << " requires total copy number " << _rtallydata.copynumvec[i]
                           << "during restart but current state only has "<<copyNum
                           <<" molecules.. Check chemistry input file. Exiting"
                           << endl;
                throw std::logic_error("Restart unsuccessful!");
            }
        }

    }
};
#endif
