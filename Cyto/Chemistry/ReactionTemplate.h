//
//  ReactionTemplate.h
//  Cyto
//
//  Created by James Komianos on 9/22/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__ReactionTemplate__
#define __Cyto__ReactionTemplate__

#include "CCylinder.h"
#include <iostream>

///Enumeration for species types
enum SpeciesType {
    BULK, DIFFUSING, FILAMENT, BOUND, PLUSEND, MINUSEND
};
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
    
    FilamentReactionDirection _direction; ///< direction of this reaction (forward, backward, in place)
    
public:
    ///Default constructor and destructor
    ReactionFilamentTemplate() {}
    ~ReactionFilamentTemplate() {}

    ///Add this chemical reaction to a CCylinder. Adds all extension and retraction callbacks needed
    virtual void addReaction(CCylinder* cc) = 0;
    
    ///Add this chemical reaction cross two CCylinders.
    ///@note assumes cc1 and cc2 are in order, that is, cc2 is the next cylinder after cc1 
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2 ) = 0;
    
};

class PolymerizationTemplate : public ReactionFilamentTemplate {
    
public:
    ///Default constructor and destructor
    PolymerizationTemplate() {}
    ~PolymerizationTemplate() {}
    
    ///to be implemented
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};

class DepolymerizationTemplate : public ReactionFilamentTemplate {
    
public:
    ///Default constructor and destructor
    DepolymerizationTemplate() {}
    ~DepolymerizationTemplate() {}
    
    ///to be implemented
    virtual void addReaction(CCylinder* cc);
    virtual void addReaction(CCylinder* cc1, CCylinder* cc2);
};


///ReactionBulkTemplate is a class to store bulk chemical reaction information read from an input file
class ReactionBulkTemplate;


#endif /* defined(__Cyto__ReactionTemplate__) */
