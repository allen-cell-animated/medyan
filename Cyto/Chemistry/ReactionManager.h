//
//  ReactionTemplate.h
//  Cyto
//
//  Created by James Komianos on 9/22/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__ReactionTemplate__
#define __Cyto__ReactionTemplate__

#include <iostream>
#include <vector>

#include "common.h"
#include "NeighborListContainer.h"
#include "CompartmentContainer.h"

///Enumeration for direction of reaction
enum FilamentReactionDirection {
    FORWARD, BACKWARD, INPLACE
};

///FORWARD DECLARATIONS
class SubSystem;
class Filament;
class CCylinder;

///InternalFilamentRxnManager is a class to store filament chemical reaction information read from an input file
/*!
 *  InternalFilamentRxnManager is used to store a filament reaction. It contains vectors of tuples that represent
 *  the position in the CMonomer in which the species is stored (for products and reactants), as well as the rate
 *  of the reaction The integer value that is the position of the species in the CMonomer vector is held by the ChemManager.
 *
 *  @note if the species is a bulk or diffusing species, the integer molecule value in the tuple stored in the SpeciesNamesDB.
 *
 *  This class also has functions to add the filament reaction to a CCylinder, as well as add a connection
 *  reaction between two neighboring CCylinders
 */

class InternalFilamentRxnManager {
    
    friend class SimpleManagerImpl;
    
protected:
    static SubSystem* _ps;
    
    ///Species identifier vectors
    vector<tuple<int,SpeciesType>> _reactants; ///< reactants in this reaction
    vector<tuple<int,SpeciesType>> _products; ///< products in this reaction
    
    float _rate; ///< rate of reaction
    
public:
    ///Default constructor and destructor
    InternalFilamentRxnManager(vector<tuple<int, SpeciesType>> reactants,
                               vector<tuple<int, SpeciesType>> products,
                               float rate) : _reactants(reactants), _products(products), _rate(rate) {}
    ~InternalFilamentRxnManager() {}

    ///Add this chemical reaction to a CCylinder. Adds all extension and retraction callbacks needed
    virtual void addReaction(CCylinder* cc) = 0;
    
    ///Add this chemical reaction cross two CCylinders.
    ///@note assumes cc1 and cc2 are in order, that is, cc2 is the next cylinder after cc1 
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) = 0;
    
    ///return compartment grid key
    CompartmentGridKey compartmentGridKey() {return CompartmentGridKey();}
};

///Template for polymerization at plus end
class PolyPlusEndManager : public InternalFilamentRxnManager {

public:
    ///Default constructor and destructor
    PolyPlusEndManager(vector<tuple<int, SpeciesType>> reactants,
                       vector<tuple<int, SpeciesType>> products,
                       float rate) : InternalFilamentRxnManager(reactants, products, rate) {}
    ~PolyPlusEndManager() {}
    
    ///to be implemented
    virtual void addReaction(CCylinder* cc);
    
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};

///Template for polymerization at minus end
class PolyMinusEndManager : public InternalFilamentRxnManager {
    
public:
    ///Default constructor and destructor
    PolyMinusEndManager(vector<tuple<int, SpeciesType>> reactants,
                        vector<tuple<int, SpeciesType>> products,
                        float rate) : InternalFilamentRxnManager(reactants, products, rate) {}
    ~PolyMinusEndManager() {}
    
    ///to be implemented
    virtual void addReaction(CCylinder* cc);
    
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};


///Template for depolymerization at plus end
class DepolyPlusEndManager : public InternalFilamentRxnManager {
    
public:
    ///Default constructor and destructor
    DepolyPlusEndManager(vector<tuple<int, SpeciesType>> reactants,
                         vector<tuple<int, SpeciesType>> products,
                         float rate) : InternalFilamentRxnManager(reactants, products, rate) {}
    ~DepolyPlusEndManager() {}
    
    ///to be implemented
    virtual void addReaction(CCylinder* cc);
    
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};

///Template for depolymerization at minus end
class DepolyMinusEndManager : public InternalFilamentRxnManager {
    
public:
    ///Default constructor and destructor
    DepolyMinusEndManager(vector<tuple<int, SpeciesType>> reactants,
                          vector<tuple<int, SpeciesType>> products,
                          float rate) : InternalFilamentRxnManager(reactants, products, rate) {}
    ~DepolyMinusEndManager() {}
    
    ///to be implemented
    virtual void addReaction(CCylinder* cc);
    
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};

class BasicBindingManager : public InternalFilamentRxnManager {
    
public:
    ///default constructor and destructor
    BasicBindingManager(vector<tuple<int, SpeciesType>> reactants,
                        vector<tuple<int, SpeciesType>> products,
                        float rate) : InternalFilamentRxnManager(reactants, products, rate) {}
    ~BasicBindingManager() {}
    
    virtual void addReaction(CCylinder* cc);
    
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) {};
};

class UnbindingManager : public InternalFilamentRxnManager {
    
public:
    ///default constructor and destructor
    UnbindingManager(vector<tuple<int, SpeciesType>> reactants,
                     vector<tuple<int, SpeciesType>> products,
                     float rate) : InternalFilamentRxnManager(reactants, products, rate) {}
    ~UnbindingManager() {}
    
    virtual void addReaction(CCylinder* cc);
    
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) {};

};

class MotorWalkFManager : public InternalFilamentRxnManager {
    
public:
    ///default constructor and destructor
    MotorWalkFManager(vector<tuple<int, SpeciesType>> reactants,
                      vector<tuple<int, SpeciesType>> products,
                      float rate) : InternalFilamentRxnManager(reactants, products, rate) {}
    ~MotorWalkFManager() {}
    
    virtual void addReaction(CCylinder* cc);
    
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);

};

class MotorWalkBManager : public InternalFilamentRxnManager {
    
public:
    ///default constructor and destructor
    MotorWalkBManager(vector<tuple<int, SpeciesType>> reactants,
                      vector<tuple<int, SpeciesType>> products,
                      float rate) : InternalFilamentRxnManager(reactants, products, rate) {}
    ~MotorWalkBManager() {}
    
    virtual void addReaction(CCylinder* cc);
    
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
    
};


/// CrossFilamentRxnManager is a class to store cross-filament reactions, including linker and motor binding

/*!
 *  CrossFilamentRxnManager is used to store a cross-filament reaction. It contains vectors of tuples that represent
 *  the position in the CMonomer in which the species is stored (for products and reactants), as well as the rate
 *  of the reaction and direction. The integer value that is the position of the species in the CMonomer vector
 *  is held by the ChemManager. Also contains the range of this reaction.
 *
 *  Also a subclass of CylinderNLContainer, contains a cylinder neighbors list of active reactions
 *
 *  @note if the species is a bulk or diffusing species, the integer molecule value in the tuple stored in the SpeciesNamesDB.
 *
 *  This class also has functions to add the cross filament reaction to two CCylinders
 */

class CrossFilamentRxnManager : public CylinderNLContainer {

    friend class SimpleManagerImpl;
    
protected:
    static SubSystem* _ps;
    
    ///Species identifier vectors
    vector<tuple<int,SpeciesType>> _reactants; ///< reactants in this reaction
    vector<tuple<int,SpeciesType>> _products; ///< products in this reaction
    
    float _rate; ///< rate of reaction
    float _rMin; ///< minimum reaction range
    float _rMax; ///< maximum reaction range
    
    int _reactionID; ///Unique ID of this reaction

public:
    ///Default constructor and destructor
    CrossFilamentRxnManager(vector<tuple<int, SpeciesType>> reactants,
                            vector<tuple<int, SpeciesType>> products,
                            float rate, float rMin, float rMax, int ID)
                            : CylinderNLContainer(rMax, rMin, true),
                              _reactants(reactants), _products(products), _rate(rate), _rMin(rMin), _rMax(rMax), _reactionID(ID) {}
    ~CrossFilamentRxnManager() {}
    
    ///add this chemical reaction to two ccylinders if within the reaction range
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) = 0;
    
    ///return compartment grid key
    CompartmentGridKey compartmentGridKey() {return CompartmentGridKey();}
    
    ///Getters for rmin and rmax
    float getRMin() {return _rMin;}
    float getRMax() {return _rMax;}
    
    ///Get reactionID
    int getReactionID() {return _reactionID;}
    
};

class LinkerBindingManager : public CrossFilamentRxnManager {
    
public:
    ///default constructor and destructor
    LinkerBindingManager(vector<tuple<int, SpeciesType>> reactants,
                    vector<tuple<int, SpeciesType>> products,
                    float rate, float rMin, float rMax, int ID)
                    : CrossFilamentRxnManager(reactants, products, rate, rMin, rMax, ID) {}
    ~LinkerBindingManager() {}
    
    ///to be implemented
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};

class MotorBindingManager : public CrossFilamentRxnManager {
    
public:
    ///default constructor and destructor
    MotorBindingManager(vector<tuple<int, SpeciesType>> reactants,
                         vector<tuple<int, SpeciesType>> products,
                         float rate, float rMin, float rMax, int ID)
                         : CrossFilamentRxnManager(reactants, products, rate, rMin, rMax, ID) {}
    ~MotorBindingManager() {}
    
    ///to be implemented
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};



#endif /* defined(__Cyto__ReactionTemplate__) */
