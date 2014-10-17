//
//  ReactionTemplate.h
//  Cyto
//
//  Created by James Komianos on 9/22/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__ReactionTemplate__
#define __Cyto__ReactionTemplate__

#include "Filament.h"

#include <iostream>

///Enumeration for direction of reaction
enum FilamentReactionDirection {
    FORWARD, BACKWARD, INPLACE
};

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
    
protected:
    ///Species identifier vectors
    std::vector<std::tuple<int,SpeciesType>> _reactants; ///< reactants in this reaction
    std::vector<std::tuple<int,SpeciesType>> _products; ///< products in this reaction
    
    float _rate; ///< rate of reaction
    
public:
    ///Default constructor and destructor
    ReactionFilamentTemplate(std::vector<std::tuple<int, SpeciesType>> reactants,
                             std::vector<std::tuple<int, SpeciesType>> products,
                             float rate) : _reactants(reactants), _products(products), _rate(rate) {}
    ~ReactionFilamentTemplate() {}

    ///Add this chemical reaction to a CCylinder. Adds all extension and retraction callbacks needed
    virtual void addReaction(CCylinder* cc, Filament* pf) = 0;
    
    ///Add this chemical reaction cross two CCylinders.
    ///@note assumes cc1 and cc2 are in order, that is, cc2 is the next cylinder after cc1 
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2, Filament* pf ) = 0;
    
};

///Template for polymerization at plus end
class PolymerizationPlusEndTemplate : public ReactionFilamentTemplate {

public:
    ///Default constructor and destructor
    PolymerizationPlusEndTemplate(std::vector<std::tuple<int, SpeciesType>> reactants,
                                  std::vector<std::tuple<int, SpeciesType>> products,
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
    PolymerizationMinusEndTemplate(std::vector<std::tuple<int, SpeciesType>> reactants,
                                   std::vector<std::tuple<int, SpeciesType>> products,
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
    DepolymerizationPlusEndTemplate(std::vector<std::tuple<int, SpeciesType>> reactants,
                                    std::vector<std::tuple<int, SpeciesType>> products,
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
    DepolymerizationMinusEndTemplate(std::vector<std::tuple<int, SpeciesType>> reactants,
                                     std::vector<std::tuple<int, SpeciesType>> products,
                                     float rate) : ReactionFilamentTemplate(reactants, products, rate) {}
    ~DepolymerizationMinusEndTemplate() {}
    
    ///to be implemented
    virtual void addReaction(CCylinder* cc, Filament* pf);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2, Filament* pf);
};

class BindingTemplate : public ReactionFilamentTemplate {
    
public:
    ///default constructor and destructor
    BindingTemplate(std::vector<std::tuple<int, SpeciesType>> reactants,
                            std::vector<std::tuple<int, SpeciesType>> products,
                            float rate) : ReactionFilamentTemplate(reactants, products, rate) {}
    ~BindingTemplate() {}
    
    virtual void addReaction(CCylinder* cc1, Filament* pf);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2, Filament* pf) {};
};

class UnbindingTemplate : public ReactionFilamentTemplate {
    
public:
    ///default constructor and destructor
    UnbindingTemplate(std::vector<std::tuple<int, SpeciesType>> reactants,
                           std::vector<std::tuple<int, SpeciesType>> products,
                           float rate) : ReactionFilamentTemplate(reactants, products, rate) {}
    ~UnbindingTemplate() {}
    
    virtual void addReaction(CCylinder* cc1, Filament* pf);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2, Filament* pf) {};

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
    
protected:
    ///Species identifier vectors
    std::vector<std::tuple<int,SpeciesType>> _reactants; ///< reactants in this reaction
    std::vector<std::tuple<int,SpeciesType>> _products; ///< products in this reaction
    
    float _rate; ///< rate of reaction
    float _rMin; ///< minimum reaction range
    float _rMax; ///< maximum reaction range

public:
    ///Default constructor and destructor
    ReactionCrossFilamentTemplate(std::vector<std::tuple<int, SpeciesType>> reactants,
                             std::vector<std::tuple<int, SpeciesType>> products,
                             float rate, float rMin, float rMax) : _reactants(reactants), _products(products), _rate(rate), _rMin(rMin), _rMax(rMax) {}
    ~ReactionCrossFilamentTemplate() {}
    
    ///add this chemical reaction to two ccylinders if within the reaction range
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2) = 0;
};

class LinkerBindingTemplate : public ReactionCrossFilamentTemplate {
    
public:
    ///default constructor and destructor
    LinkerBindingTemplate(std::vector<std::tuple<int, SpeciesType>> reactants,
                    std::vector<std::tuple<int, SpeciesType>> products,
                    float rate, float rMin, float rMax) : ReactionCrossFilamentTemplate(reactants, products, rate, rMin, rMax) {}
    ~LinkerBindingTemplate() {}
    
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};

class MotorBindingTemplate : public ReactionCrossFilamentTemplate {
    
public:
    ///default constructor and destructor
    MotorBindingTemplate(std::vector<std::tuple<int, SpeciesType>> reactants,
                          std::vector<std::tuple<int, SpeciesType>> products,
                          float rate, float rMin, float rMax) : ReactionCrossFilamentTemplate(reactants, products, rate, rMin, rMax) {}
    ~MotorBindingTemplate() {}
    
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};


///ReactionBulkTemplate is a class to store bulk chemical reaction information read from an input file
class ReactionBulkTemplate;


#endif /* defined(__Cyto__ReactionTemplate__) */
