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
#include "CompartmentContainer.h"

///Enumeration for direction of reaction
enum FilamentReactionDirection {
    FORWARD, BACKWARD, INPLACE
};

///FORWARD DECLARATIONS
class SubSystem;
class Filament;
class CCylinder;

///ReactionFilamentTemplate is a class to store filament chemical reaction information read from an input file
/*!
 *  ReactionFilamentTemplate is used to store a filament reaction. It contains vectors of tuples that represent
 *  the position in the CMonomer in which the species is stored (for products and reactants), as well as the rate
 *  of the reaction and direction. The integer value that is the position of the species in the CMonomer vector
 *  is held by the chemical initializer.
 *
 * @note if the species is a bulk or diffusing species, the integer molecule value in the tuple stored in the SpeciesNamesDB.
 *
 *  This class also has functions to add the filament reaction to a CCylinder, as well as add a connection
 *  reaction between two neighboring CCylinders
 */

class ReactionFilamentTemplate {
    
    friend class SimpleManagerImpl;
    
protected:
    static SubSystem* _ps;
    
    ///Species identifier vectors
    vector<tuple<int,SpeciesType>> _reactants; ///< reactants in this reaction
    vector<tuple<int,SpeciesType>> _products; ///< products in this reaction
    
    float _rate; ///< rate of reaction
    
public:
    ///Default constructor and destructor
    ReactionFilamentTemplate(vector<tuple<int, SpeciesType>> reactants,
                             vector<tuple<int, SpeciesType>> products,
                             float rate) : _reactants(reactants), _products(products), _rate(rate) {}
    ~ReactionFilamentTemplate() {}

    ///Add this chemical reaction to a CCylinder. Adds all extension and retraction callbacks needed
    virtual void addReaction(CCylinder* cc, Filament* pf) = 0;
    
    ///Add this chemical reaction cross two CCylinders.
    ///@note assumes cc1 and cc2 are in order, that is, cc2 is the next cylinder after cc1 
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2, Filament* pf ) = 0;
    
    ///return compartment grid key
    CompartmentGridKey compartmentGridKey() {return CompartmentGridKey();}
};

///Template for polymerization at plus end
class PolymerizationPlusEndTemplate : public ReactionFilamentTemplate {

public:
    ///Default constructor and destructor
    PolymerizationPlusEndTemplate(vector<tuple<int, SpeciesType>> reactants,
                                  vector<tuple<int, SpeciesType>> products,
                                  float rate) : ReactionFilamentTemplate(reactants, products, rate) {}
    ~PolymerizationPlusEndTemplate() {}
    
    ///to be implemented
    virtual void addReaction(CCylinder* cc, Filament* pf);
    
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2, Filament* pf);
};

///Template for polymerization at minus end
class PolymerizationMinusEndTemplate : public ReactionFilamentTemplate {
    
public:
    ///Default constructor and destructor
    PolymerizationMinusEndTemplate(vector<tuple<int, SpeciesType>> reactants,
                                   vector<tuple<int, SpeciesType>> products,
                                   float rate) : ReactionFilamentTemplate(reactants, products, rate) {}
    ~PolymerizationMinusEndTemplate() {}
    
    ///to be implemented
    virtual void addReaction(CCylinder* cc, Filament* pf);
    
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2, Filament* pf);
};


///Template for depolymerization at plus end
class DepolymerizationPlusEndTemplate : public ReactionFilamentTemplate {
    
public:
    ///Default constructor and destructor
    DepolymerizationPlusEndTemplate(vector<tuple<int, SpeciesType>> reactants,
                                    vector<tuple<int, SpeciesType>> products,
                                    float rate) : ReactionFilamentTemplate(reactants, products, rate) {}
    ~DepolymerizationPlusEndTemplate() {}
    
    ///to be implemented
    virtual void addReaction(CCylinder* cc, Filament* pf);
    
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2, Filament* pf);
};

///Template for depolymerization at minus end
class DepolymerizationMinusEndTemplate : public ReactionFilamentTemplate {
    
public:
    ///Default constructor and destructor
    DepolymerizationMinusEndTemplate(vector<tuple<int, SpeciesType>> reactants,
                                     vector<tuple<int, SpeciesType>> products,
                                     float rate) : ReactionFilamentTemplate(reactants, products, rate) {}
    ~DepolymerizationMinusEndTemplate() {}
    
    ///to be implemented
    virtual void addReaction(CCylinder* cc, Filament* pf);
    
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2, Filament* pf);
};

class BasicBindingTemplate : public ReactionFilamentTemplate {
    
public:
    ///default constructor and destructor
    BasicBindingTemplate(vector<tuple<int, SpeciesType>> reactants,
                         vector<tuple<int, SpeciesType>> products,
                         float rate) : ReactionFilamentTemplate(reactants, products, rate) {}
    ~BasicBindingTemplate() {}
    
    virtual void addReaction(CCylinder* cc1, Filament* pf);
    
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2, Filament* pf) {};
};

class UnbindingTemplate : public ReactionFilamentTemplate {
    
public:
    ///default constructor and destructor
    UnbindingTemplate(vector<tuple<int, SpeciesType>> reactants,
                           vector<tuple<int, SpeciesType>> products,
                           float rate) : ReactionFilamentTemplate(reactants, products, rate) {}
    ~UnbindingTemplate() {}
    
    virtual void addReaction(CCylinder* cc1, Filament* pf);
    
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2, Filament* pf) {};

};

class MotorWalkingForwardTemplate : public ReactionFilamentTemplate {
    
public:
    ///default constructor and destructor
    MotorWalkingForwardTemplate(vector<tuple<int, SpeciesType>> reactants,
                                vector<tuple<int, SpeciesType>> products,
                                float rate) : ReactionFilamentTemplate(reactants, products, rate) {}
    ~MotorWalkingForwardTemplate() {}
    
    virtual void addReaction(CCylinder* cc1, Filament* pf);
    
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2, Filament* pf);

};

class MotorWalkingBackwardTemplate : public ReactionFilamentTemplate {
    
public:
    ///default constructor and destructor
    MotorWalkingBackwardTemplate(vector<tuple<int, SpeciesType>> reactants,
                                vector<tuple<int, SpeciesType>> products,
                                float rate) : ReactionFilamentTemplate(reactants, products, rate) {}
    ~MotorWalkingBackwardTemplate() {}
    
    virtual void addReaction(CCylinder* cc1, Filament* pf);
    
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2, Filament* pf);
    
};


///ReactionCrossFilamentTemplate is a class to store cross-filament reactions, including linker and motor binding/unbinding

/*!
 *  ReactionCrossFilamentTemplate is used to store a cross-filament reaction. It contains vectors of tuples that represent
 *  the position in the CMonomer in which the species is stored (for products and reactants), as well as the rate
 *  of the reaction and direction. The integer value that is the position of the species in the CMonomer vector
 *  is held by the chemical initializer. Also contains the range of this reaction.
 *
 * @note if the species is a bulk or diffusing species, the integer molecule value in the tuple stored in the SpeciesNamesDB.
 *
 *  This class also has functions to add the filament reaction to two CCylinders
 */

class ReactionCrossFilamentTemplate {

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
    ReactionCrossFilamentTemplate(vector<tuple<int, SpeciesType>> reactants,
                             vector<tuple<int, SpeciesType>> products, float rate, float rMin, float rMax, int ID)
                             : _reactants(reactants), _products(products), _rate(rate), _rMin(rMin), _rMax(rMax), _reactionID(ID) {}
    ~ReactionCrossFilamentTemplate() {}
    
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

class LinkerBindingTemplate : public ReactionCrossFilamentTemplate {
    
public:
    ///default constructor and destructor
    LinkerBindingTemplate(vector<tuple<int, SpeciesType>> reactants,
                    vector<tuple<int, SpeciesType>> products,
                    float rate, float rMin, float rMax, int ID)
                    : ReactionCrossFilamentTemplate(reactants, products, rate, rMin, rMax, ID) {}
    ~LinkerBindingTemplate() {}
    
    ///to be implemented
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};

class MotorBindingTemplate : public ReactionCrossFilamentTemplate {
    
public:
    ///default constructor and destructor
    MotorBindingTemplate(vector<tuple<int, SpeciesType>> reactants,
                         vector<tuple<int, SpeciesType>> products,
                         float rate, float rMin, float rMax, int ID)
                         : ReactionCrossFilamentTemplate(reactants, products, rate, rMin, rMax, ID) {}
    ~MotorBindingTemplate() {}
    
    ///to be implemented
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};



#endif /* defined(__Cyto__ReactionTemplate__) */
