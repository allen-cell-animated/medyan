
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

#include "CBranchingPoint.h"

#include "ChemCallbacks.h"
#include "CompartmentGrid.h"
#include "CCylinder.h"
#include "CMonomer.h"

CBranchingPoint::CBranchingPoint(short branchType, Compartment* c,
                                 CCylinder* cc1, CCylinder* cc2, int position)

    : CBound(cc1->getType(), c, cc1, cc2, position, 0), _branchType(branchType) {

    //Find species on cylinder that should be marked
    SpeciesBound* sb1 = _cc1->getCMonomer(_position1)->speciesBrancher(branchType);
    SpeciesBound* se1 = _cc1->getCMonomer(_position1)->speciesBound(
                        SysParams::Chemistry().brancherBoundIndex[_filamentType]);

//    //@{
//    SpeciesBound* BL1 = _cc1->getCMonomer(_position1)->speciesBound(
//            SysParams::Chemistry().linkerBoundIndex[_filamentType]);
//    SpeciesBound* BL2 = _cc2->getCMonomer(_position2)->speciesBound(
//            SysParams::Chemistry().linkerBoundIndex[_filamentType]);
//    SpeciesBound* BB1 = _cc1->getCMonomer(_position1)->speciesBound(
//            SysParams::Chemistry().brancherBoundIndex[_filamentType]);
//    SpeciesBound* BB2 = _cc2->getCMonomer(_position2)->speciesBound(
//            SysParams::Chemistry().brancherBoundIndex[_filamentType]);
//    SpeciesBound* BM1 = _cc1->getCMonomer(_position1)->speciesBound(
//            SysParams::Chemistry().motorBoundIndex[_filamentType]);
//    SpeciesBound* BM2 = _cc2->getCMonomer(_position2)->speciesBound(
//            SysParams::Chemistry().motorBoundIndex[_filamentType]);
//    SpeciesBound* sm1 = _cc1->getCMonomer(_position1)->speciesMotor(0);
//    SpeciesBound* sm2 = _cc2->getCMonomer(_position2)->speciesMotor(0);
//    SpeciesBound* sl1 = _cc1->getCMonomer(_position1)->speciesLinker(0);
//    SpeciesBound* sl2 = _cc2->getCMonomer(_position2)->speciesLinker(0);
//
//    std::cout<<"BranchingPoint "<<cc1->getCylinder()->getID()<<" "<<_position1<<" "<<cc2->getCylinder()->getID()<<" "<<
//             ""<<_position2<<" branchType "<<branchType<<endl;
//    std::cout<<"Motor "<<sm1->getN()<<" "<<sm2->getN()<<" BOUND "<<BM1->getN()<<" "<<BM2->getN()<<endl;
//    std::cout<<"Linker "<<sl1->getN()<<" "<<sl2->getN()<<" BOUND "<<BL1->getN()<<" "<<BL2->getN()<<endl;
//    std::cout<<"Brancher "<<sb1->getN()<<" BOUND "<<BB1->getN()<<" "<<BB2->getN()<<endl;
//    //@}


    //mark species
    assert(areEqual(sb1->getN(), 0.0) && areEqual(se1->getN(), 1.0) &&
           "Major bug: Brancher binding to an occupied site.");
        
    sb1->up(); se1->down();

    //attach this branchpoint to the species
    setFirstSpecies(sb1);
}

CBranchingPoint::~CBranchingPoint() {
    
    //remove the unbinding reaction
    _cc1->removeInternalReaction(_offRxn);
    
}

void CBranchingPoint::createOffReaction(ReactionBase* onRxn, SubSystem* ps){
    
    //first, find the correct diffusing or bulk species
    RSpecies** rs = onRxn->rspecies();
    Species* sfb = &(rs[SPECIESB_BINDING_INDEX]->getSpecies());
    
    //create the reaction species
    CMonomer* m = _cc1->getCMonomer(_position1);
    vector<Species*> os = {m->speciesBrancher(_branchType),
                           m->speciesBound(SysParams::Chemistry().brancherBoundIndex[_filamentType]), sfb};
    
    //create reaction, add to cylinder
    ReactionBase* offRxn =
    new Reaction<BUNBINDINGREACTANTS,BUNBINDINGPRODUCTS>(os, _offRate);
    
    offRxn->setReactionType(ReactionType::BRANCHUNBINDING);
    
    //add the unbinding reaction and callback
    BranchingPointUnbindingCallback bcallback(_pBranchingPoint, ps);
    ConnectionBlock rcb(offRxn->connect(bcallback,false));
    
    setOffReaction(offRxn);
    _cc1->addInternalReaction(offRxn);
    
}
