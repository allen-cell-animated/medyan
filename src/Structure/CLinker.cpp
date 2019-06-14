
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

#include "CLinker.h"

#include "ChemCallbacks.h"
#include "CCylinder.h"
#include "Cylinder.h"

CLinker::CLinker(short linkerType, Compartment* c,
                 CCylinder* cc1, CCylinder* cc2, int position1, int position2)

    : CBound(cc1->getType(), c, cc1, cc2, position1, position2) {
        
    //Find species on cylinder that should be marked
    SpeciesBound* sl1 = _cc1->getCMonomer(_position1)->speciesLinker(linkerType);
    SpeciesBound* sl2 = _cc2->getCMonomer(_position2)->speciesLinker(linkerType);

    SpeciesBound* se1 = _cc1->getCMonomer(_position1)->speciesBound(
                        SysParams::Chemistry().linkerBoundIndex[_filamentType]);
    SpeciesBound* se2 = _cc2->getCMonomer(_position2)->speciesBound(
                        SysParams::Chemistry().linkerBoundIndex[_filamentType]);
////    //@{
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
//    SpeciesBound* sb1 = _cc1->getCMonomer(_position1)->speciesBrancher(0);
//    SpeciesBound* sb2 = _cc2->getCMonomer(_position2)->speciesBrancher(0);
//    std::cout<<"Linker "<<cc1->getCylinder()->getID()<<" "<<_position1<<" "<<cc2->getCylinder()->getID()<<" "<<
//            ""<<_position2<<" linkerType "<<linkerType<<endl;
//	cout<<"Linker cIndices "<<cc1->getCylinder()->_dcIndex<<" "<<cc2->getCylinder()
//			->_dcIndex<<endl;
//    std::cout<<"Species Bound "<<se1->getN()<<" "<<se2->getN()<<endl;
   /* std::cout<<"Motor "<<sm1->getN()<<" "<<sm2->getN()<<" BOUND "<<BM1->getN()<<" "
                                                                                 ""<<BM2->getN()<<endl;*/
//    std::cout<<"Linker "<<sl1->getN()<<" "<<sl2->getN()<<" BOUND "<<BL1->getN()<<" "<<BL2->getN()<<endl;
    /*std::cout<<"Brancher "<<sb1->getN()<<" "<<sb2->getN()<<" BOUND "<<BB1->getN()<<"
    "<<BB2->getN()<<endl;*/
//    std::cout<<sl1->getN()<<" "<<sl2->getN()<<" "<<se1->getN()<<" "<<se2->getN()<<endl;

    /*for(auto c:Cylinder::getCylinders()){
        std::cout<<c->getID()<<" "<<c->getMCylinder()->getLength()<<" ";
    }
    std::cout<<endl;*/
//    //@}

#ifdef DETAILEDOUTPUT
    std::cout<<"Chosen sites Cyl1 "<<cc1->getCylinder()->getID()<<" bs1 "<<_position1<<" "
            "Cyl2 "<<cc2->getCylinder()->getID()<<" bs2 "<<_position2<<endl;
#endif

    //mark species
        
    assert(areEqual(sl1->getN(), (floatingpoint)0.0) && areEqual(sl2->getN(), (floatingpoint)0.0) &&
           areEqual(se1->getN(), (floatingpoint)1.0) && areEqual(se2->getN(), (floatingpoint)1.0) &&
           "Major bug: Linker binding to an occupied site.");
    
    sl1->up(); sl2->up();
    se1->down(); se2->down();
    
    //attach this linker to the species
    setFirstSpecies(sl1);
    setSecondSpecies(sl2);
}

CLinker::~CLinker() {

    //remove the off reaction
    _cc1->removeCrossCylinderReaction(_cc2, _offRxn);

}

void CLinker::createOffReaction(ReactionBase* onRxn, SubSystem* ps) {
    
    RSpecies** rs = onRxn->rspecies();
    vector<Species*> os;
    
    //copy into offspecies vector
    os.push_back(_firstSpecies);
    os.push_back(_secondSpecies);
    
    os.push_back(&rs[SPECIESL_BINDING_INDEX]->getSpecies());
    
    Species* empty1 = _cc1->getCMonomer(_position1)->speciesBound(
                      SysParams::Chemistry().linkerBoundIndex[_filamentType]);
    Species* empty2 = _cc2->getCMonomer(_position2)->speciesBound(
                      SysParams::Chemistry().linkerBoundIndex[_filamentType]);
    
    os.push_back(empty1);
    os.push_back(empty2);
    
    ReactionBase* offRxn =
    new Reaction<LMUNBINDINGREACTANTS,LMUNBINDINGPRODUCTS>(os, _offRate);
    offRxn->setReactionType(ReactionType::LINKERUNBINDING);
    
    //set gnum of offreaction
    // Dissipation
    if(SysParams::Chemistry().dissTracking){

    floatingpoint gnum = onRxn->getGNumber();
    offRxn->setGNumber(-gnum);

    //set hrcdid of offreaction
    string hrcdid = onRxn->getHRCDID();
    offRxn->setHRCDID(hrcdid + "off");
    }

    //Attach the callback to the off reaction, add it
    LinkerUnbindingCallback lcallback(_pLinker, ps);
    ConnectionBlock rcb(offRxn->connect(lcallback,false));
    
    _cc1->addCrossCylinderReaction(_cc2, offRxn);
    setOffReaction(offRxn);
}
