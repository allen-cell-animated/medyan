
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

#include "CMotorGhost.h"

#include "ChemCallbacks.h"
#include "CCylinder.h"

CMotorGhost::CMotorGhost(short motorType, Compartment* c,
                         CCylinder* cc1, CCylinder* cc2, int position1, int position2)

    : CBound(cc1->getType(), c, cc1, cc2, position1, position2) {
    
    //Find species on cylinder that should be marked
    SpeciesBound* sm1 = _cc1->getCMonomer(_position1)->speciesMotor(motorType);
    SpeciesBound* sm2 = _cc2->getCMonomer(_position2)->speciesMotor(motorType);
    SpeciesBound* se1 = _cc1->getCMonomer(_position1)->speciesBound(
                        SysParams::Chemistry().motorBoundIndex[_filamentType]);
    SpeciesBound* se2 = _cc2->getCMonomer(_position2)->speciesBound(
                        SysParams::Chemistry().motorBoundIndex[_filamentType]);
//    //@{
    /*SpeciesBound* BL1 = _cc1->getCMonomer(_position1)->speciesBound(
            SysParams::Chemistry().linkerBoundIndex[_filamentType]);
    SpeciesBound* BL2 = _cc2->getCMonomer(_position2)->speciesBound(
            SysParams::Chemistry().linkerBoundIndex[_filamentType]);
    SpeciesBound* BB1 = _cc1->getCMonomer(_position1)->speciesBound(
            SysParams::Chemistry().brancherBoundIndex[_filamentType]);
    SpeciesBound* BB2 = _cc2->getCMonomer(_position2)->speciesBound(
            SysParams::Chemistry().brancherBoundIndex[_filamentType]);*/
//    SpeciesBound* BM1 = _cc1->getCMonomer(_position1)->speciesBound(
//            SysParams::Chemistry().motorBoundIndex[_filamentType]);
//    SpeciesBound* BM2 = _cc2->getCMonomer(_position2)->speciesBound(
//            SysParams::Chemistry().motorBoundIndex[_filamentType]);
    /*SpeciesBound* sl1 = _cc1->getCMonomer(_position1)->speciesLinker(0);
    SpeciesBound* sl2 = _cc2->getCMonomer(_position2)->speciesLinker(0);
    SpeciesBound* sb1 = _cc1->getCMonomer(_position1)->speciesBrancher(0);
    SpeciesBound* sb2 = _cc2->getCMonomer(_position2)->speciesBrancher(0);*/
//    std::cout<<"Motor "<<cc1->getCylinder()->getID()<<" "<<_position1<<" "<<cc2->getCylinder()->getID()<<" "<<
//             ""<<_position2<<" MotorType "<<motorType<< endl;
//    cout<<"Motor cIndices "<<cc1->getCylinder()->_dcIndex<<" "<<cc2->getCylinder()
//    ->_dcIndex<<endl;
//        std::cout<<"Motor "<<sm1->getN()<<" "<<sm2->getN()<<" BOUND "<<BM1->getN()<<" "<<BM2->getN()<<endl;
//        std::cout<<"Linker "<<sl1->getN()<<" "<<sl2->getN()<<" BOUND "<<BL1->getN()<<" "<<BL2->getN()<<endl;
//        std::cout<<"Brancher "<<sb1->getN()<<" "<<sb2->getN()<<" BOUND "<<BB1->getN()<<" "<<BB2->getN()<<endl;
//    for(auto c:Cylinder::getCylinders()){
//        std::cout<<c->getId()<<" "<<c->getMCylinder()->getLength()<<" ";
//    }
//    std::cout<<endl;
//    //@}
#ifdef DETAILEDOUTPUT
    std::cout<<"Chosen sites Cyl1 "<<cc1->getCylinder()->getId()<<" bs1 "<<_position1<<" "
            "Cyl2 "<<cc2->getCylinder()->getId()<<" bs2 "<<_position2<<endl;
#endif
    /*std::cout<<"Motor binding "<<cc1->getCylinder()->getID()<<" "<<_position1<<" "
             <<cc2->getCylinder()->getID()<<" "<<_position2<<" MotorType "<<motorType<< endl;
    cout<<"species statuses (bound) "<<se1->getN()<<" "<<se2->getN()<<" (motor) "<<sm1->getN()<<" "
									""<<sm2->getN()<<endl;
	auto c1coord = cc1->getCompartment()->coordinates();
	auto c2coord = cc2->getCompartment()->coordinates();
    cout<<"coordinate Cmp "<<c1coord[0]<<" "<<c1coord[1]<<" "<<c1coord[2]<<" "
																		   ""<<c2coord[0]<<" "<<c2coord[1]<<" "<<c2coord[2]<<endl;*/
    //mark species
    assert(areEqual(sm1->getN(), 0.0) && areEqual(sm2->getN(), 0.0) &&
           areEqual(se1->getN(), 1.0) && areEqual(se2->getN(), 1.0) &&
           "Major bug: Motor binding to an occupied site.");
        
    sm1->up(); sm2->up();
    se1->down();
    se2->down();
        
    //attach this motor to the species
    setFirstSpecies(sm1);
    setSecondSpecies(sm2);
}

CMotorGhost::~CMotorGhost() {
    #ifdef CHECKRXN
	cout<<"Remove cross mID destructor "<<_pMotorGhost->getId()<<" Cylinders "
														""<<_cc1->getCylinder()
			->getId()<<" "<<_cc2->getCylinder()->getId()<<" offRxn "<<_offRxn<<endl;
    #endif
    
    //remove the off reaction
    _cc1->removeCrossCylinderReaction(_cc2, _offRxn);
    
}

void CMotorGhost::createOffReaction(ReactionBase* onRxn, SubSystem* ps) {
    
    RSpecies** rs = onRxn->rspecies();
    vector<Species*> os;
    
    //copy into offspecies vector
    os.push_back(_firstSpecies);
    os.push_back(_secondSpecies);
    
    os.push_back(&rs[SPECIESM_BINDING_INDEX]->getSpecies());
    
    Species* empty1 = _cc1->getCMonomer(_position1)->speciesBound(
                      SysParams::Chemistry().motorBoundIndex[_filamentType]);
    Species* empty2 = _cc2->getCMonomer(_position2)->speciesBound(
                      SysParams::Chemistry().motorBoundIndex[_filamentType]);
    
    os.push_back(empty1);
    os.push_back(empty2);
    
    ReactionBase* offRxn =
    new Reaction<LMUNBINDINGREACTANTS,LMUNBINDINGPRODUCTS>(os, _offRate);

    offRxn->setReactionType(ReactionType::MOTORUNBINDING);
    // Dissipation
    if(SysParams::Chemistry().dissTracking){
    floatingpoint gnum = onRxn->getGNumber();
    offRxn->setGNumber(-gnum);
    
    //set hrcdid of offreaction
    string hrcdid = onRxn->getHRCDID();
    offRxn->setHRCDID(hrcdid + "off");
    }
    //Attach the callback to the off reaction, add it
    MotorUnbindingCallback mcallback(_pMotorGhost, ps);
    ConnectionBlock rcb(offRxn->connect(mcallback,false));
	#ifdef CHECKRXN
	cout<<"Add cross mID "<<_pMotorGhost->getId()<<" Cylinders "<<_cc1->getCylinder()
			->getId()<<" "<<_cc2->getCylinder()->getId()<<" offRxn "<<_offRxn<<endl;
	#endif
    _cc1->addCrossCylinderReaction(_cc2, offRxn);
    setOffReaction(offRxn);
    #ifdef CHECKRXN
    cout<<"create MG off rxn "<<offRxn<<"with RNode "<<offRxn->getRnode()<<" with mID "
    <<_pMotorGhost->getId()<<endl;
	#endif
}

void CMotorGhost::moveMotorHead(CCylinder* cc,
                                short oldPosition,
                                short newPosition,
                                short motorType,
                                short boundType,
                                SubSystem* ps) {
    
    auto smOld = cc->getCMonomer(oldPosition)->speciesMotor(motorType);
    auto smNew = cc->getCMonomer(newPosition)->speciesMotor(motorType);

    auto seNew = cc->getCMonomer(newPosition)->speciesBound(boundType);
    
    ReactionBase* newOffRxn;
    
    // Dissipation
    string hrcdid = "NA";
    floatingpoint gnum = 0.0;
    if(SysParams::Chemistry().dissTracking){
        hrcdid = _offRxn->getHRCDID();
        gnum = _offRxn->getGNumber();
    }
    
    if(getFirstSpecies() == smOld) {
        
        _position1 = newPosition;
        
        setFirstSpecies(smNew);
        
        //change off reaction to include new species
        Species* smOther = _secondSpecies;
        Species* seOther = _cc2->getCMonomer(_position2)->speciesBound(
                           SysParams::Chemistry().motorBoundIndex[_filamentType]);
        
        Species* sbd = &(_offRxn->rspecies()[SPECIESM_UNBINDING_INDEX]->getSpecies());
        
        //create new reaction
        newOffRxn = new Reaction<LMUNBINDINGREACTANTS, LMUNBINDINGPRODUCTS>
                          ({smNew, smOther, sbd, seNew, seOther}, _offRate);
    }
    else {
        _position2 = newPosition;
        
        setSecondSpecies(smNew);
        
        //change off reaction to include new species
        Species* smOther = _firstSpecies;
        Species* seOther = _cc1->getCMonomer(_position1)->speciesBound(
                           SysParams::Chemistry().motorBoundIndex[_filamentType]);
        
        Species* sbd = &(_offRxn->rspecies()[SPECIESM_UNBINDING_INDEX]->getSpecies());
        
        //create new reaction
        newOffRxn = new Reaction<LMUNBINDINGREACTANTS, LMUNBINDINGPRODUCTS>
                          ({smOther, smNew, sbd, seNew, seOther}, _offRate);
    }
    //set new reaction type
    newOffRxn->setReactionType(ReactionType::MOTORUNBINDING);
    
    // set new reaction gnum
    
    // Dissipation
    if(SysParams::Chemistry().dissTracking){
    newOffRxn->setGNumber(gnum);
    
    //set hrcdid of offreaction
    
    newOffRxn->setHRCDID(hrcdid);
    }
    
    //attach signal
    MotorUnbindingCallback mcallback(_pMotorGhost, ps);
    ConnectionBlock rcb(newOffRxn->connect(mcallback,false));
    #ifdef CHECKRXN
	cout<<"Remove cross mID "<<_pMotorGhost->getId()<<" Cylinders "<<_cc1->getCylinder()
			->getId()<<" "<<_cc2->getCylinder()->getId()<<" offRxn "<<_offRxn<<endl;
    #endif
    //remove old reaction, add new one
    _cc1->removeCrossCylinderReaction(_cc2, _offRxn);
	#ifdef CHECKRXN
	cout<<"Add cross mID "<<_pMotorGhost->getId()<<" Cylinders "<<_cc1->getCylinder()
			->getId()<<" "<<_cc2->getCylinder()->getId()<<" offRxn "<<newOffRxn<<endl;
	#endif
    _cc1->addCrossCylinderReaction(_cc2, newOffRxn);
    
    //set new unbinding reaction
    setOffReaction(newOffRxn);
    
}

void CMotorGhost::moveMotorHead(CCylinder* oldCC,
                                CCylinder* newCC,
                                short oldPosition,
                                short newPosition,
                                short motorType,
                                short boundType,
                                SubSystem* ps) {

    
    auto smOld = oldCC->getCMonomer(oldPosition)->speciesMotor(motorType);
    auto smNew = newCC->getCMonomer(newPosition)->speciesMotor(motorType);
    
    auto seNew = newCC->getCMonomer(newPosition)->speciesBound(boundType);
    
    ReactionBase* newOffRxn;
    
    
    // Dissipation
    string hrcdid = "NA";
    floatingpoint gnum = 0.0;
    if(SysParams::Chemistry().dissTracking){
    hrcdid = _offRxn->getHRCDID();
    gnum = _offRxn->getGNumber();
    }
    if(getFirstSpecies() == smOld) {
        
        _position1 = newPosition;
        
        setFirstSpecies(smNew);
        
        //change off reaction to include new species
        Species* smOther = _secondSpecies;
        Species* seOther = _cc2->getCMonomer(_position2)->speciesBound(
                           SysParams::Chemistry().motorBoundIndex[_filamentType]);
        
        Species* sbd = &(_offRxn->rspecies()[SPECIESM_UNBINDING_INDEX]->getSpecies());
        
        //create new reaction
        newOffRxn = new Reaction<LMUNBINDINGREACTANTS, LMUNBINDINGPRODUCTS>
                          ({smNew, smOther, sbd, seNew, seOther}, _offRate);
        
        //remove old off reaction
        #ifdef CHECKRXN
	    cout<<"Remove cross mID "<<_pMotorGhost->getId()<<" Cylinders "<<_cc1->getCylinder()
			    ->getId()<<" "<<_cc2->getCylinder()->getId()<<" offRxn "<<_offRxn<<endl;
        #endif
        _cc1->removeCrossCylinderReaction(_cc2, _offRxn);
        
        //set new ccylinders
        setFirstCCylinder(newCC);
    }
    else {
        _position2 = newPosition;
        
        setSecondSpecies(smNew);
        
        //change off reaction to include new species
        Species* smOther = _firstSpecies;
        Species* seOther = _cc1->getCMonomer(_position1)->speciesBound(
                           SysParams::Chemistry().motorBoundIndex[_filamentType]);
        
        Species* sbd = &(_offRxn->rspecies()[SPECIESM_UNBINDING_INDEX]->getSpecies());
        
        //create new reaction
        newOffRxn = new Reaction<LMUNBINDINGREACTANTS, LMUNBINDINGPRODUCTS>
                         ({smOther, smNew, sbd, seNew, seOther}, _offRate);

        #ifdef CHECKRXN
        cout<<"Remove cross mID "<<_pMotorGhost->getId()<<" Cylinders "<<_cc1->getCylinder()
                ->getId()<<" "<<_cc2->getCylinder()->getId()<<" offRxn "<<_offRxn<<endl;
		#endif
        //remove old off reaction
        _cc1->removeCrossCylinderReaction(_cc2, _offRxn);
        
        //set new ccylinders
        setSecondCCylinder(newCC);
    }
    
    //set new reaction type
    newOffRxn->setReactionType(ReactionType::MOTORUNBINDING);
    
    // Dissipation
    if(SysParams::Chemistry().dissTracking){
    // set new reaction gnum
    newOffRxn->setGNumber(gnum);
    
    //set hrcdid of offreaction
    
    newOffRxn->setHRCDID(hrcdid);
    }
    //attach signal
    MotorUnbindingCallback mcallback(_pMotorGhost, ps);
    
    ConnectionBlock rcb(newOffRxn->connect(mcallback,false));
	#ifdef CHECKRXN
	cout<<"Add cross mID "<<_pMotorGhost->getId()<<" Cylinders "<<_cc1->getCylinder()
			->getId()<<" "<<_cc2->getCylinder()->getId()<<" offRxn "<<newOffRxn<<endl;
	#endif
    //add new
    _cc1->addCrossCylinderReaction(_cc2, newOffRxn);
    
    //set new unbinding reaction
    setOffReaction(newOffRxn);
    
}


void CMotorGhost::printReaction(){
    cout<<_offRxn->getHRCDID()<<endl;
}
