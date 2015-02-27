
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

#ifndef M3SYM_RateChanger_h
#define M3SYM_RateChanger_h

#include "common.h"

/// Used to change Filament reaction rates based on forces in the network
/*!
 *  The FilamentRateChanger class is an abstract class which allows
 *  for Filament rate changing based on a given force. Different
 *  implementations of this class will have different rate changing models, 
 *  and will all implement the changeRate() function.
 */

class FilamentRateChanger {
    
public:
    /// Change the reaction rate based on a bare rate and given force.
    virtual float changeRate(float bareRate, double force) = 0;
};

/// Used to change Linker reaction rates based on forces in the network
/*!
 *  The LinkerRateChanger class is an abstract class which allows
 *  for Linker rate changing based on a given force. Different 
 *  implementations of this class will have different rate changing models, 
 *  and will all implement the changeRate() function.
 */
class LinkerRateChanger {
    
protected:
    short _linkerType; ///< This linker type
    
public:
    LinkerRateChanger(short linkerType) : _linkerType(linkerType) {}
    
    /// Change the reaction rate based on a bare rate and given force.
    virtual float changeRate(float bareRate, double force) = 0;
};

/// Used to change MotorGhost reaction rates based on forces in the network
/*!
 *  The MotorRateChanger class is an abstract class which allows
 *  for MotorGhost rate changing based on a given force. Different 
 *  implementations of this class will have different rate changing models, 
 *  and will all implement the changeRate() function.
 */
class MotorRateChanger {
    
protected:
    short _motorType; ///< This motor type
    
public:
    MotorRateChanger(short motorType) : _motorType(motorType) {}
    
    /// Change the reaction rate based on a bare rate,
    /// number of heads, and given force.
    virtual float changeRate(float bareRate, int numHeads, double force) = 0;
};


#endif
