
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

#include <cmath>

#include "MotorGhost.h"

#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"
#include "ChemRNode.h"

#include "GController.h"
#include "SysParams.h"
#include "MathFunctions.h"
#include "Mechanics/CUDAcommon.h"
#include "Rand.h"

using namespace mathfunc;

void MotorGhost::updateCoordinate() {
    
    auto x1 = _c1->getFirstBead()->vcoordinate();
    auto x2 = _c1->getSecondBead()->vcoordinate();
    auto x3 = _c2->getFirstBead()->vcoordinate();
    auto x4 = _c2->getSecondBead()->vcoordinate();
    
    auto m1 = midPointCoordinate(x1, x2, _position1);
    auto m2 = midPointCoordinate(x3, x4, _position2);
    
    coordinate = midPointCoordinate(m1, m2, 0.5);
}


MotorGhost::MotorGhost(Cylinder* c1, Cylinder* c2, short motorType,
                       floatingpoint position1, floatingpoint position2,
                       floatingpoint onRate, floatingpoint offRate)

    : Trackable(true, true),
      _c1(c1), _c2(c2),
      _position1(position1), _position2(position2),
      _motorType(motorType), _birthTime(tau()),
      _onRate(onRate), _offRate(offRate) {
          
    //find compartment
    updateCoordinate();
    
    try {_compartment = GController::getCompartment(coordinate);}
    catch (exception& e) {
        cout << e.what();
        
        printSelf();
        
        exit(EXIT_FAILURE);
    }
    short filType = c1->getType();
          
    int pos1 = int(position1 * SysParams::Geometry().cylinderNumMon[filType]);
    int pos2 = int(position2 * SysParams::Geometry().cylinderNumMon[filType]);
          
    //set number of heads by picking random int between maxheads and minheads
    _numHeads = Rand::randInteger(SysParams::Chemistry().motorNumHeadsMin[_motorType],
                                  SysParams::Chemistry().motorNumHeadsMax[_motorType]);
    
    if(!_unbindingChangers.empty())
        _numBoundHeads = _unbindingChangers[_motorType]->numBoundHeads(_onRate, _offRate, 0, _numHeads);
    else
        _numBoundHeads = _numHeads;
    
#ifdef CHEMISTRY
    _cMotorGhost = unique_ptr<CMotorGhost>(
    new CMotorGhost(motorType, _compartment, _c1->getCCylinder(), _c2->getCCylinder(), pos1, pos2));
    _cMotorGhost->setMotorGhost(this);
#endif

#ifdef MECHANICS
    auto x1 = _c1->getFirstBead()->vcoordinate();
    auto x2 = _c1->getSecondBead()->vcoordinate();
    auto x3 = _c2->getFirstBead()->vcoordinate();
    auto x4 = _c2->getSecondBead()->vcoordinate();
#ifdef PLOSFEEDBACK
    _mMotorGhost = unique_ptr<MMotorGhost>(
    new MMotorGhost(motorType, _numHeads, position1, position2, x1, x2, x3, x4));
    _mMotorGhost->setMotorGhost(this);
#else
    //Using _numBoundHeads vs _numHeads as the argument to the constructor here - Cal edit to produce realistic dissipation values
    _mMotorGhost = unique_ptr<MMotorGhost>(
            new MMotorGhost(motorType, _numBoundHeads, position1, position2, x1, x2, x3, x4));
    _mMotorGhost->setMotorGhost(this);
#endif
#endif
#ifdef CROSSCHECK_MOTOR
	auto b1 = _c1->getFirstBead();
	auto b2 = _c1->getSecondBead();
	auto b3 = _c2->getFirstBead();
	auto b4 = _c2->getSecondBead();

	cout<<"motor binding  mID "<<getId()<<" cIDs "<<_c1->getId()<<" "<<_c2->getId()<<" "
	                                                                                 ""<<_position1<<" "
	    <<_position2<<" cindex "<<_c1->getStableIndex()<<" "<<_c2->getStableIndex()<<" bID "<<b1->getId
			()<<" "<<b2->getId()<<" bindex "<<b1->getStableIndex()<<" "<<b2->getStableIndex()<<" "
	    <<b3->getStableIndex()<<" "<<b4->getStableIndex()<<endl;
#endif
    
}

///@note - record lifetime data here
MotorGhost::~MotorGhost() noexcept {
	#ifdef CROSSCHECK_MOTOR
	auto b1 = _c1->getFirstBead();
	auto b2 = _c1->getSecondBead();
	auto b3 = _c2->getFirstBead();
	auto b4 = _c2->getSecondBead();

	cout<<"motor un-binding  mID "<<getId()<<" cIDs "<<_c1->getId()<<" "<<_c2->getId()<<" "
	<<_position1<<" "<<_position2<<" cindex "<<_c1->getStableIndex()<<" "
	<<_c2->getStableIndex()<<" bID "<<b1->getId()<<" "<<b2->getId()<<" bindex "
	<<b1->getStableIndex()<<" "<<b2->getStableIndex()<<" "
	<<b3->getStableIndex()<<" "<<b4->getStableIndex()<<endl;
	#endif

//    floatingpoint lifetime = tau() - _birthTime;
//    
//    if(_lifetimes->getMax() > lifetime &&
//       _lifetimes->getMin() < lifetime)
//        _lifetimes->addValue(lifetime);
//        
//    
//    if(_walkLengths->getMax() > _walkLength &&
//       _walkLengths->getMin() < _walkLength)
//        _walkLengths->addValue(_walkLength);
    
}

void MotorGhost::updatePosition() {
#ifdef CROSSCHECK_IDX
    cout<<"MG position begin"<<endl;
#endif
#ifdef CHEMISTRY
    //update ccylinders
    _cMotorGhost->setFirstCCylinder(_c1->getCCylinder());
    _cMotorGhost->setSecondCCylinder(_c2->getCCylinder());
    
#endif
    //check if in same compartment
    updateCoordinate();
    Compartment* c;
    
    try {c = GController::getCompartment(coordinate);}
    catch (exception& e) {
        cout << e.what();
        
        printSelf();
        
        exit(EXIT_FAILURE);
    }
    
    if(c != _compartment) {
	    #ifdef CROSSCHECK_MOTOR
    	auto oldcoord = _compartment->coordinates();
	    auto newcoord = c->coordinates();

    	cout<<"MotorGhost  move Cmp mID "<<this->getId()<<" from Cmp ID "
    	<<_compartment->getId()<<" with Coord "<<oldcoord[0]<<" " <<oldcoord[1]<<" "
    	<<oldcoord[2]<<" to Cmp ID "
    	<<c->getId()<<" with Coord "<<newcoord[0]<<" " <<newcoord[1]<<" "
	    <<newcoord[2]<<endl;
		#endif

        mins = chrono::high_resolution_clock::now();
        
        _compartment = c;
#ifdef CHEMISTRY
        SpeciesBound* firstSpecies = _cMotorGhost->getFirstSpecies();
        SpeciesBound* secondSpecies = _cMotorGhost->getSecondSpecies();
#ifdef CHECKRXN
        cout<<"Clone _cMotorGhost with offRxn "<<_cMotorGhost->getOffReaction()<<" with "
																				 "Id "
					 ""<<_cMotorGhost->getOffReaction()->getId()<<endl;
#endif
        CMotorGhost* clone = _cMotorGhost->clone(c);
        setCMotorGhost(clone);
        
        _cMotorGhost->setFirstSpecies(firstSpecies);
        _cMotorGhost->setSecondSpecies(secondSpecies);
#endif
        mine = chrono::high_resolution_clock::now();
        chrono::duration<floatingpoint> compartment_update(mine - mins);
        CUDAcommon::tmin.timemotorupdate += compartment_update.count();
        CUDAcommon::tmin.callsmotorupdate++;
    }
    
#ifdef MECHANICS
    auto x1 = _c1->getFirstBead()->vcoordinate();
    auto x2 = _c1->getSecondBead()->vcoordinate();
    auto x3 = _c2->getFirstBead()->vcoordinate();
    auto x4 = _c2->getSecondBead()->vcoordinate();
    
    auto m1 = midPointCoordinate(x1, x2, _position1);
    auto m2 = midPointCoordinate(x3, x4, _position2);
    
    _mMotorGhost->setLength(twoPointDistance(m1, m2));
    
    //update the spring constant, based on numboundheads
    //current force
    floatingpoint force = max((floatingpoint)0.0, _mMotorGhost->stretchForce);
    
    //update number of bound heads
    if(!_unbindingChangers.empty())
        _numBoundHeads = _unbindingChangers[_motorType]->numBoundHeads(_onRate, _offRate, force, _numHeads);
    else
        _numBoundHeads = _numHeads;
    
#ifdef PLOSFEEDBACK
    _mMotorGhost->setStretchingConstant(_motorType, _numHeads);
#else
    _mMotorGhost->setStretchingConstant(_motorType, _numBoundHeads);
#endif

#endif
#ifdef CROSSCHECK_IDX
	cout<<"MG position end"<<endl;
#endif
}

/// @note - This function updates forward walking rates using the
/// stetching force in the opposite direction of the motor walk.
/// Does not consider negative forces in this direction.

/// Updates unbinding rates based on the stretch force. Does not
/// consider compression forces, only stretching.
void MotorGhost::updateReactionRates() {

#ifdef CROSSCHECK_IDX
    auto b1 = _c1->getFirstBead();
    auto b2 = _c1->getSecondBead();
    auto b3 = _c2->getFirstBead();
    auto b4 = _c2->getSecondBead();
    cout<<"MG mID "<<getId()<<" "<<_c1->getId()<<" "<<_c2->getId()<<" "<<_position1<<" "
        <<_position2<<" cindex "<<_c1->getStableIndex()<<" "<<_c2->getStableIndex()<<" bID "<<b1->getId
            ()<<" "<<b2->getId()<<" bindex "<<b1->getStableIndex()<<" "<<b2->getStableIndex()<<" "
            <<b3->getStableIndex()<<" "<<b4->getStableIndex()<<endl;
#endif

    //current force
    floatingpoint force = max<floatingpoint>((floatingpoint)0.0,
            _mMotorGhost->stretchForce);
    
    //update number of bound heads
    if(!_unbindingChangers.empty())
        _numBoundHeads = _unbindingChangers[_motorType]->numBoundHeads(_onRate, _offRate, force, _numHeads);
    else
        _numBoundHeads = _numHeads;
    
    //walking rate changer
    if(!_walkingChangers.empty()) {
        auto x1 = _c1->getFirstBead()->vcoordinate();
        auto x2 = _c1->getSecondBead()->vcoordinate();
        auto x3 = _c2->getFirstBead()->vcoordinate();
        auto x4 = _c2->getSecondBead()->vcoordinate();

	    const auto& cylinderInfoData = Cylinder::getDbData().value;

	    auto c1struct = cylinderInfoData[_c1->getStableIndex()];
	    auto c2struct = cylinderInfoData[_c2->getStableIndex()];
	    auto fType1 = c1struct.type;
	    auto fType2 = c2struct.type;

	    bool consider_passivation = false;
	    bool isc1leftofc2 = false;
	    /* isc1leftofc2 tells which of the two cylinders is close to minus end. c1 is
	     * connected to leg1 and c2 is connected to leg2 of motor.
	     * If isc1leftofc2= true, in MotorWalkingForward reactions, the reaction
	     * involving leg1 should be passivated.
	     * If isc1leftofc2 = false, in MotorWalkingBackward reactions, the reaction
	     * involving leg2 should be passivated.*/

	    if((c1struct.filamentId == c2struct.filamentId)) {
	    	auto c1posonFil = c1struct.positionOnFilament;
		    auto c2posonFil = c2struct.positionOnFilament;
            consider_passivation = abs(c1posonFil - c2posonFil) <= SysParams::Mechanics().sameFilBindSkip + 1;
		    //A distance of 3 or lesser between two cylinders on the same filament is not
		    // acceptable.
		    isc1leftofc2 = c1posonFil < c2posonFil;
	    }
/*#ifdef MOTORBIASCHECK
	    consider_passivation = false;
#endif*/

        auto mp1 = midPointCoordinate(x1, x2, _position1);
        auto mp2 = midPointCoordinate(x3, x4, _position2);

        //get component of force in direction of forward walk for C1, C2
        vector<floatingpoint> motorC1Direction =
        twoPointDirection(mp1, mp2);
        
        vector<floatingpoint> motorC2Direction =
        twoPointDirection(mp2, mp1);
        
        vector<floatingpoint> c1Direction = twoPointDirection(x2,x1);
        vector<floatingpoint> c2Direction = twoPointDirection(x4,x3);
        
        floatingpoint forceDotDirectionC1 = force * dotProduct(motorC1Direction, c1Direction);
        floatingpoint forceDotDirectionC2 = force * dotProduct(motorC2Direction, c2Direction);
        
        //WALKING REACTIONS
        Species* s1 = _cMotorGhost->getFirstSpecies();
        Species* s2 = _cMotorGhost->getSecondSpecies();
        
        for(auto r : s1->getRSpecies().reactantReactions()) {
            
            if(r->getReactionType() == ReactionType::MOTORWALKINGFORWARD) {

            	if(consider_passivation && isc1leftofc2){
		            #ifdef CROSSCHECK_MOTOR
		            auto b1 = _c1->getFirstBead();
		            auto b2 = _c1->getSecondBead();
		            auto b3 = _c2->getFirstBead();
		            auto b4 = _c2->getSecondBead();
		            cout<<"motor-walking  c1 mID "<<getId()<<" cIDs "<<_c1->getId()<<" "
		                                                                             ""<<_c2->getId()<<" "
		                                                                                               ""<<_position1<<" "
		                <<_position2<<" cindex "<<_c1->getStableIndex()<<" "<<_c2->getStableIndex()<<" bID "<<b1->getId
				            ()<<" "<<b2->getId()<<" bindex "<<b1->getStableIndex()<<" "<<b2->getStableIndex()<<" "
		                <<b3->getStableIndex()<<" "<<b4->getStableIndex()<<endl;
		            #endif

		            float c1lastbindingsite = float(*(SysParams::Chemistry()
				            .bindingSites[fType1].end() -1))/float(SysParams::Geometry()
				            		.cylinderNumMon[fType1]);
		            #ifdef CROSSCHECK_MOTOR
		            cout<<"c1lastbindingsite "<<c1lastbindingsite<<" "<<
		            *(SysParams::Chemistry().bindingSites[fType1].begin())<<" "
		            << *(SysParams::Chemistry().bindingSites[fType1].end()-1)<<" "<<
		            SysParams::Geometry().cylinderNumMon[fType1]<<endl;
					#endif
		            if(areEqual(c1lastbindingsite,_position1))
			            r->setRateMulFactor(0.0f, ReactionBase::MOTORWALKCONSTRAINTFACTOR);
            	}
            	else{
		            r->setRateMulFactor(1.0f, ReactionBase::MOTORWALKCONSTRAINTFACTOR);
		            float newRate =
				            _walkingChangers[_motorType]->
						            changeRate(_cMotorGhost->getOnRate(),
						                       _cMotorGhost->getOffRate(),
						                       _numHeads, max<floatingpoint>((floatingpoint)0.0, forceDotDirectionC1));
		            if(SysParams::RUNSTATE==false){
			            newRate=0.0;}
#ifdef DETAILEDOUTPUT
		            std::cout<<"Motor WF1 f "<<force<<" Rate "<<newRate<<" "<<coordinate[0]<<" "
                        ""<<coordinate[1]<<" "<<coordinate[2]<<" Fdirn "<<
                         forceDotDirectionC2<<" NH "<<_numHeads<<endl;
#endif
		            r->setBareRate(newRate);
	            }/*if(consider_passivation && isc1leftofc2)*/
	            r->updatePropensity();
            }/*MOTORWALKINGFORWARD*/
            else if(r->getReactionType() == ReactionType::MOTORWALKINGBACKWARD) {

	            if(consider_passivation && !isc1leftofc2) {

		            auto c1firstbindingsite = float(*(SysParams::Chemistry()
				            .bindingSites[fType1].begin()))/
		                                      float(SysParams::Geometry()
		                                      .cylinderNumMon[fType1]);
		            if(areEqual(c1firstbindingsite,_position1))
			            r->setRateMulFactor(0.0f, ReactionBase::MOTORWALKCONSTRAINTFACTOR);

	            }else{
		            r->setRateMulFactor(1.0f, ReactionBase::MOTORWALKCONSTRAINTFACTOR);
		            float newRate =
				            _walkingChangers[_motorType]->
						            changeRate(_cMotorGhost->getOnRate(),
						                       _cMotorGhost->getOffRate(),
						                       _numHeads, max<floatingpoint>((floatingpoint)0.0, -forceDotDirectionC1));

		            if(SysParams::RUNSTATE==false){ newRate=0.0;}
#ifdef DETAILEDOUTPUT
		            std::cout<<"Motor WB1 f "<<force<<" Rate "<<newRate<<" "<<coordinate[0]<<" "
                        ""<<coordinate[1]<<" "<<coordinate[2]<<" Fdirn "<<
                         forceDotDirectionC2<<" NH "<<_numHeads<<endl;
#endif
		            r->setBareRate(newRate);
	            }
	            r->updatePropensity();
            }
        }
        for(auto r : s2->getRSpecies().reactantReactions()) {
            
            if(r->getReactionType() == ReactionType::MOTORWALKINGFORWARD) {

	            if(consider_passivation && !isc1leftofc2) {
		            #ifdef CROSSCHECK_MOTOR
		            auto b1 = _c1->getFirstBead();
		            auto b2 = _c1->getSecondBead();
		            auto b3 = _c2->getFirstBead();
		            auto b4 = _c2->getSecondBead();
		            cout<<"motor-walking  c2 mID "<<getId()<<" cIDs "<<_c1->getId()<<" "
		                                                                             ""<<_c2->getId()<<" "
		                                                                                               ""<<_position1<<" "
		                <<_position2<<" cindex "<<_c1->getStableIndex()<<" "<<_c2->getStableIndex()<<" bID "<<b1->getId
				            ()<<" "<<b2->getId()<<" bindex "<<b1->getStableIndex()<<" "<<b2->getStableIndex()<<" "
		                <<b3->getStableIndex()<<" "<<b4->getStableIndex()<<endl;
		            #endif

		            float c2lastbindingsite = float(*(SysParams::Chemistry()
				            .bindingSites[fType2].end()-1))
				            		/float(SysParams::Geometry().cylinderNumMon[fType2]);
		            #ifdef CROSSCHECK_MOTOR
		            cout<<"c2lastbindingsite "<<c2lastbindingsite<<" "<<
		                *(SysParams::Chemistry().bindingSites[fType1].begin())<<" "
		                << *(SysParams::Chemistry().bindingSites[fType1].end()-1)<<" "<<
		                SysParams::Geometry().cylinderNumMon[fType1]<<endl;
		            #endif
		            if(areEqual(c2lastbindingsite,_position2))
			            r->setRateMulFactor(0.0f, ReactionBase::MOTORWALKCONSTRAINTFACTOR);

	            }else {
		            r->setRateMulFactor(1.0f, ReactionBase::MOTORWALKCONSTRAINTFACTOR);
		            if(r->isPassivated()) {r->activateReaction();}
		            float newRate =
				            _walkingChangers[_motorType]->
						            changeRate(_cMotorGhost->getOnRate(),
						                       _cMotorGhost->getOffRate(),
						                       _numHeads,
						                       max<floatingpoint>((floatingpoint) 0.0,
						                                          forceDotDirectionC2));
		            if (SysParams::RUNSTATE == false) { newRate = 0.0; }

#ifdef DETAILEDOUTPUT
					std::cout<<"Motor WF2 f "<<force<<" Rate "<<newRate<<" "<<coordinate[0]<<" "
							""<<coordinate[1]<<" "<<coordinate[2]<<" Fdirn "<<
							forceDotDirectionC2<<" NH "<<_numHeads<<endl;
#endif
		            r->setBareRate(newRate);
	            }
	            r->updatePropensity();
            }/*MOTORWALKINGFORWARD*/
            else if(r->getReactionType() == ReactionType::MOTORWALKINGBACKWARD) {

	            if(consider_passivation && isc1leftofc2) {
		            auto c2firstbindingsite = float(*(SysParams::Chemistry()
				            .bindingSites[fType2].begin()))/float(SysParams::Geometry()
				                                      .cylinderNumMon[fType2]);
		            if(areEqual(c2firstbindingsite,_position2))
			            r->setRateMulFactor(0.0f, ReactionBase::MOTORWALKCONSTRAINTFACTOR);
	            }
	            else{
		            r->setRateMulFactor(1.0f, ReactionBase::MOTORWALKCONSTRAINTFACTOR);
	            	if(r->isPassivated()) {r->activateReaction();}
		            float newRate =
				            _walkingChangers[_motorType]->
						            changeRate(_cMotorGhost->getOnRate(),
						                       _cMotorGhost->getOffRate(),
						                       _numHeads, max<floatingpoint>((floatingpoint)0.0, -forceDotDirectionC2));
		            if(SysParams::RUNSTATE==false)
		            { newRate=0.0;}
#ifdef DETAILEDOUTPUT
		            std::cout<<"Motor WB2 f "<<force<<" Rate "<<newRate<<" "<<coordinate[0]<<" "
                        ""<<coordinate[1]<<" "<<coordinate[2]<<" Fdirn "<<
                         forceDotDirectionC2<<" NH "<<_numHeads<<endl;
#endif
		            r->setBareRate(newRate);
	            }
	            r->updatePropensity();
            }
        }
    }
    
    //unbinding rate changer
    if(!_unbindingChangers.empty()) {

        //get the unbinding reaction
        ReactionBase* offRxn = _cMotorGhost->getOffReaction();

        #ifdef CHECKRXN
        cout<<"update MG off rxn "<<offRxn<<" with RNode "<<offRxn->getRnode()<<" with "
																				"mID "<<this->getId()<<endl;
        #endif
        
        //change the rate
        float newRate =
        _unbindingChangers[_motorType]->
        changeRate(_cMotorGhost->getOnRate(), _cMotorGhost->getOffRate(), _numHeads, force);
#ifdef DETAILEDOUTPUT
        std::cout<<"Motor UB f "<<force<<" Rate "<<newRate<<" "<<coordinate[0]<<" "
                ""<<coordinate[1]<<" "
                ""<<coordinate[2]<<endl;
#endif
        offRxn->setBareRate(newRate);
        offRxn->activateReaction();
    }
#ifdef CROSSCHECK
    cout<<"MG mID "<<getId()<<" "<<_c1->getId()<<" "<<_c2->getId()<<" "<<_position1<<" "
        <<_position2<<" cindex "<<_c1->getStableIndex()<<" "<<_c2->getStableIndex()<<" bID "<<b1->getId
            ()<<" "<<b2->getId()<<" bindex "<<b1->getStableIndex()<<" "<<b2->getStableIndex()<<" "
        <<b3->getStableIndex()<<" "<<b4->getStableIndex()<<endl;
#endif
}

void MotorGhost::moveMotorHead(Cylinder* c,
                               floatingpoint oldPosition, floatingpoint newPosition,
                               short boundType, SubSystem* ps) {
#ifdef CROSSCHECK_MOTOR
	//Print coordinate
	auto x1 = _c1->getFirstBead()->vcoordinate();
	auto x2 = _c1->getSecondBead()->vcoordinate();
	auto x3 = _c2->getFirstBead()->vcoordinate();
	auto x4 = _c2->getSecondBead()->vcoordinate();
	auto mp1 = midPointCoordinate(x1, x2, _position1);
	auto mp2 = midPointCoordinate(x3, x4, _position2);
	auto motorcoord = midPointCoordinate(mp1,mp2,0.5);

    auto b1 = _c1->getFirstBead();
    auto b2 = _c1->getSecondBead();
    auto b3 = _c2->getFirstBead();
    auto b4 = _c2->getSecondBead();


	cout<<"motor-walking  mID "<<getId()<<" cIDs "<<_c1->getId()<<" "<<_c2->getId()<<" "
                                                                              ""<<_position1<<" "
        <<_position2<<" cindex "<<_c1->getStableIndex()<<" "<<_c2->getStableIndex()<<" bID "<<b1->getId
            ()<<" "<<b2->getId()<<" bindex "<<b1->getStableIndex()<<" "<<b2->getStableIndex()<<" "
        <<b3->getStableIndex()<<" "<<b4->getStableIndex()<<" move "<<c->getId()<<" old "<<oldPosition<<" new "
        <<newPosition<<endl;
	cout<<"mcoord before "<<mp1[0]<<" "<<mp1[1]<<" "<<mp1[2]<<" "<<mp2[0]<<" "<<mp2[1]<<" "
																			  ""<<mp2[2]<<endl;

    if(c->getId() == _c1->getId() && newPosition == _position1)
    	cout<<"c1 walking on itself"<<endl;
    else if(c->getId() == _c2->getId() && newPosition == _position2)
	    cout<<"c2 walking on itself"<<endl;
#endif

    //shift the position of one side of the motor
    floatingpoint shift =  newPosition - oldPosition;
    
    //shift the head
    if(c == _c1) {
        _position1 += shift;
    }
    else {
        _position2 += shift;
    }
    short filType = c->getType();
    
    //record walk length
    _walkLength += shift * SysParams::Geometry().cylinderSize[filType];
    
#ifdef CHEMISTRY
    short oldpos = int (oldPosition * SysParams::Geometry().cylinderNumMon[filType]);
    short newpos = int (newPosition * SysParams::Geometry().cylinderNumMon[filType]);
    
    _cMotorGhost->moveMotorHead(c->getCCylinder(), oldpos, newpos,
                                _motorType, boundType, ps);
#if CROSSCHECK_MOTOR
	x1 = _c1->getFirstBead()->vcoordinate();
	x2 = _c1->getSecondBead()->vcoordinate();
	x3 = _c2->getFirstBead()->vcoordinate();
	x4 = _c2->getSecondBead()->vcoordinate();
	mp1 = midPointCoordinate(x1, x2, _position1);
	mp2 = midPointCoordinate(x3, x4, _position2);
	auto DeltaL = twoPointDistance(mp1, mp2) - _mMotorGhost->getEqLength();

    //(If motor is stretched) stretchForce > 0 => DeltaL > 0
    //(If motor is contracted) stretchForce < 0 => DeltaL < 0
    //Neither stretchForce = 0 => DeltaL = 0
    bool stretchstate = _mMotorGhost->stretchForce > 0.0;
    bool stretchstatenow = DeltaL > 0.0;

    bool unstretched = _mMotorGhost->stretchForce == 0.0;
    bool unstretchednow = DeltaL == 0.0;

    //stretch - stretch
    if(stretchstate == true && stretchstatenow == true)
        CUDAcommon::mwalk.stretchtostretch++;
    //stretch - contract
    else if(stretchstate == true && stretchstatenow == false )
        CUDAcommon::mwalk.stretchtocontract++;
    //contract - contract
    else if(stretchstate == false && stretchstatenow == false )
        CUDAcommon::mwalk.contracttocontract++;
    //contract - stretch
    else if(stretchstate == false && stretchstatenow == true)
        CUDAcommon::mwalk.contracttostretch++;
    //equib - contract
    else if(unstretched == true && stretchstatenow == false)
        CUDAcommon::mwalk.equibtocontract++;
    //equib - stretch
    else if(unstretched == true && stretchstatenow == true)
        CUDAcommon::mwalk.equibtostretch++;
    //equib - equib
    else if(unstretched == true && unstretchednow == true)
        CUDAcommon::mwalk.equibtoequib++;

	cout<<"mcoord after "<<mp1[0]<<" "<<mp1[1]<<" "<<mp1[2]<<" "<<mp2[0]<<" "<<mp2[1]<<" "<<mp2[2]<<endl;
#endif


#endif
    
}

void MotorGhost::moveMotorHead(Cylinder* oldC, Cylinder* newC,
                               floatingpoint oldPosition, floatingpoint newPosition,
                               short boundType, SubSystem* ps) {
#ifdef CROSSCHECK_MOTOR
	//Print coordinate
	auto x1 = _c1->getFirstBead()->vcoordinate();
	auto x2 = _c1->getSecondBead()->vcoordinate();
	auto x3 = _c2->getFirstBead()->vcoordinate();
	auto x4 = _c2->getSecondBead()->vcoordinate();
	auto mp1 = midPointCoordinate(x1, x2, _position1);
	auto mp2 = midPointCoordinate(x3, x4, _position2);
	auto motorcoord = midPointCoordinate(mp1,mp2,0.5);

    auto b1 = _c1->getFirstBead();
    auto b2 = _c1->getSecondBead();
    auto b3 = _c2->getFirstBead();
    auto b4 = _c2->getSecondBead();

    cout<<"motor-walking  mID "<<getId()<<" "<<_c1->getId()<<" "<<_c2->getId()<<" "<<_position1<<" "
        <<_position2<<" cindex "<<_c1->getStableIndex()<<" "<<_c2->getStableIndex()<<" bID "<<b1->getId
            ()<<" "<<b2->getId()<<" bindex "<<b1->getStableIndex()<<" "<<b2->getStableIndex()<<" "
        <<b3->getStableIndex()<<" "<<b4->getStableIndex()<<" move oldC "<<oldC->getId()<<" newC "
        <<newC->getId() <<" old " <<oldPosition<<" new "<<newPosition<<endl;

	if(newC->getId() == _c1->getId() && newPosition == _position1)
		cout<<"c1 walking on itself"<<endl;
	else if(newC->getId() == _c2->getId() && newPosition == _position2)
		cout<<"c2 walking on itself"<<endl;
#endif
    //shift the head
    if(oldC == _c1) {
        _position1 = newPosition;
        _c1 = newC;
    }
    else {
        _position2 = newPosition;
        _c2 = newC;
    }
    short filType = _c1->getType();
    
    //record walk length
    _walkLength += (1-oldPosition + newPosition) * SysParams::Geometry().cylinderSize[filType];
    
#ifdef CHEMISTRY
    short oldpos = int (oldPosition * SysParams::Geometry().cylinderNumMon[filType]);
    short newpos = int (newPosition * SysParams::Geometry().cylinderNumMon[filType]);
    
    _cMotorGhost->moveMotorHead(oldC->getCCylinder(), newC->getCCylinder(),
                                oldpos, newpos, _motorType, boundType, ps);

#if CROSSCHECK_MOTOR
	x1 = _c1->getFirstBead()->vcoordinate();
	x2 = _c1->getSecondBead()->vcoordinate();
	x3 = _c2->getFirstBead()->vcoordinate();
	x4 = _c2->getSecondBead()->vcoordinate();
	mp1 = midPointCoordinate(x1, x2, _position1);
	mp2 = midPointCoordinate(x3, x4, _position2);
	auto DeltaL = twoPointDistance(mp1, mp2) - _mMotorGhost->getEqLength();

	//(If motor is stretched) stretchForce > 0 => DeltaL > 0
	//(If motor is contracted) stretchForce < 0 => DeltaL < 0
	//Neither stretchForce = 0 => DeltaL = 0
	bool stretchstate = _mMotorGhost->stretchForce > 0.0;
	bool stretchstatenow = DeltaL > 0.0;

	bool unstretched = _mMotorGhost->stretchForce == 0.0;
	bool unstretchednow = DeltaL == 0.0;

	//stretch - stretch
	if(stretchstate == true && stretchstatenow == true)
		CUDAcommon::mwalk.stretchtostretch++;
		//stretch - contract
	else if(stretchstate == true && stretchstatenow == false )
		CUDAcommon::mwalk.stretchtocontract++;
		//contract - contract
	else if(stretchstate == false && stretchstatenow == false )
		CUDAcommon::mwalk.contracttocontract++;
		//contract - stretch
	else if(stretchstate == false && stretchstatenow == true)
		CUDAcommon::mwalk.contracttostretch++;
		//equib - contract
	else if(unstretched == true && stretchstatenow == false)
		CUDAcommon::mwalk.equibtocontract++;
		//equib - stretch
	else if(unstretched == true && stretchstatenow == true)
		CUDAcommon::mwalk.equibtostretch++;
		//equib - equib
	else if(unstretched == true && unstretchednow == true)
		CUDAcommon::mwalk.equibtoequib++;
#endif
#endif
}


void MotorGhost::printSelf()const {
    
    cout << endl;
    
    cout << "MotorGhost: ptr = " << this << endl;
    cout << "Motor type = " << _motorType << ", Motor ID = " << getId() << endl;
    cout << "Coordinates = " << coordinate[0] << ", " << coordinate[1] << ", " << coordinate[2] << endl;
    
    cout << "Position on first cylinder (floatingpoint) = " << _position1 << endl;
    cout << "Position on second cylinder (floatingpoint) = " << _position2 << endl;
    
    cout << "Number of heads = " << _numHeads << endl;
    cout << "Birth time = " << _birthTime << endl;
    
    cout << endl;
    
#ifdef CHEMISTRY
    cout << "Associated species 1 = " << _cMotorGhost->getFirstSpecies()->getName()
    << " , copy number = " << _cMotorGhost->getFirstSpecies()->getN()
    << " , position on first cylinder (int) = " << _cMotorGhost->getFirstPosition() << endl;
    
    cout << "Associated species 2 = " << _cMotorGhost->getSecondSpecies()->getName()
    << " , copy number = " << _cMotorGhost->getSecondSpecies()->getN()
    << " , position on second cylinder (int) = " << _cMotorGhost->getSecondPosition() << endl;
#endif
    
    cout << endl;
    
    cout << "Associated cylinders (one and two): " << endl;
    _c1->printSelf();
    _c2->printSelf();
    
    cout << endl;
}

species_copy_t MotorGhost::countSpecies(const string& name) {
    
    species_copy_t copyNum = 0;
    
    for(auto m : getElements()) {
        
        auto s = m->getCMotorGhost()->getFirstSpecies();
        string sname = SpeciesNamesDB::removeUniqueFilName(s->getName());
        
        if(sname == name)
            copyNum += s->getN();
    }
    return copyNum;
}

vector<MotorRateChanger*> MotorGhost::_unbindingChangers;
vector<MotorRateChanger*> MotorGhost::_walkingChangers;

Histogram* MotorGhost::_lifetimes;
Histogram* MotorGhost::_walkLengths;

