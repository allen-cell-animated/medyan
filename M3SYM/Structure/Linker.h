
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_Linker_h
#define M3SYM_Linker_h

#include "common.h"

#include "Composite.h"
#include "CLinker.h"
#include "MLinker.h"
#include "Movable.h"
#include "Reactable.h"

//FORWARD DECLARATIONS
class Cylinder;

/// Linker class is a container to store a [MLinker] (@ref MLinker) and [CLinker](@ref CLinker)
/*!
 * Linker class is used to manage and store a [MLinker] (@ref MLinker) and [CLinker](@ref CLinker).
 * Upon intialization, both of these components are created. Extending the [Movable](@ref Movable) and [Reactable] (@ref Reactable)
 * classes, the Linker can update its position and reactions according to mechanical equilibration.
 */

class Linker : public Composite, public Movable, public Reactable {

private:
    unique_ptr<MLinker> _mLinker; ///< Pointer to MLinker
    unique_ptr<CLinker> _cLinker; ///< Pointer to CLinker
    
    Cylinder* _c1; ///< First cylinder the linker is bound to
    Cylinder* _c2; ///< Second cylinder the linker is bound to
    
    double _position1; ///< Position on first cylinder
    double _position2; ///< Position on second cylinder
    
    short _linkerType; ///< Integer specifying the type of linker
    int _linkerID; ///< Integer ID of this specific linker
    
    float _birthTime; ///Birth time of this linker
    
    Compartment* _compartment; ///< Compartment that this linker is in
    
public:
    vector<double> coordinate; ///< coordinate of midpoint, updated with updatePosition()
    
    Linker(Cylinder* c1, Cylinder* c2, short linkerType, int linkerID, double position1, double position2, bool creation);
    ~Linker();
    
    //@{
    ///Get attached cylinder
    Cylinder* getFirstCylinder() {return _c1;}
    Cylinder* getSecondCylinder() {return _c2;}
    //@}
    
    /// Set CLinker
    void setCLinker(CLinker* cLinker) {_cLinker = unique_ptr<CLinker>(cLinker);}
    /// Get CLinker
    CLinker* getCLinker() {return _cLinker.get();}
    
    /// Get MLinker
    MLinker* getMLinker() {return _mLinker.get();}
    
    //@{
    /// Linker position management
    double getFirstPosition() {return _position1;}
    void setFirstPosition(double position1) {_position1 = position1;}
    
    double getSecondPosition() {return _position2;}
    void setSecondPosition(double position2) {_position2 = position2;}
    //@}
    
    //@{
    /// Get linker parameter
    short getLinkerType() {return _linkerType;}
    int getLinkerID() {return _linkerID;}
    //@}
    
    /// Update the position of this Linker
    /// @note - changes compartment of clinker if needed
    virtual void updatePosition();
    
    /// Update the reaction rates of this linker
    virtual void updateReactionRates();
    
};


#endif 
