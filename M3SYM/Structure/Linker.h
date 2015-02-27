
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

#ifndef M3SYM_Linker_h
#define M3SYM_Linker_h

#include "common.h"

#include "LinkerDB.h"

#include "Composite.h"
#include "CLinker.h"
#include "MLinker.h"
#include "Movable.h"
#include "Reactable.h"
#include "RateChangerImpl.h"

//FORWARD DECLARATIONS
class Cylinder;
class DRController;

/// A container to store a MLinker and CLinker.
/*!
 * Linker class is used to manage and store a MLinker and CLinker. Upon intialization,
 * both of these components are created. Extending the Movable and Reactable classes, 
 * the Linker can update its position and reactions according to mechanical 
 * equilibration.
 */
class Linker : public Composite, public Movable, public Reactable {

friend class DRController;
    
private:
    unique_ptr<MLinker> _mLinker; ///< Pointer to mech linker
    unique_ptr<CLinker> _cLinker; ///< Pointer to chem linker
    
    Cylinder* _c1; ///< First cylinder the linker is bound to
    Cylinder* _c2; ///< Second cylinder the linker is bound to
    
    double _position1; ///< Position on first cylinder
    double _position2; ///< Position on second cylinder
    
    short _linkerType; ///< Integer specifying the type
    int _linkerID; ///< Integer ID of this specific linker, managed by LinkerDB
    
    float _birthTime; ///Birth time
    
    Compartment* _compartment; ///< Where this linker is
    
    ///For dynamic rate unbinding
    static vector<LinkerRateChanger*> _unbindingChangers;
    
public:
    vector<double> coordinate;
        ///< coordinate of midpoint, updated with updatePosition()
    
    Linker(Cylinder* c1, Cylinder* c2, short linkerType,
           double position1 = 0.5, double position2 = 0.5, bool creation = false);
    virtual ~Linker() noexcept;
    
    //@{
    ///Get attached cylinder
    Cylinder* getFirstCylinder() {return _c1;}
    Cylinder* getSecondCylinder() {return _c2;}
    //@}
    
    /// Set chem linker
    void setCLinker(CLinker* cLinker) {_cLinker = unique_ptr<CLinker>(cLinker);}
    /// Get chem linker
    CLinker* getCLinker() {return _cLinker.get();}
    
    /// Get mech linker
    MLinker* getMLinker() {return _mLinker.get();}
    
    //@{
    /// Position management
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
    
    /// Get the birth time
    float getBirthTime() {return _birthTime;}
    
    /// Update the position
    /// @note - changes compartment if needed
    virtual void updatePosition();
    
    /// Update the reaction rates
    virtual void updateReactionRates();
    
};


#endif 
