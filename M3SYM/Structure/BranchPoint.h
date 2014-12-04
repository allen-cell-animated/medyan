
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

#ifndef M3SYM_BranchPoint_h
#define M3SYM_BranchPoint_h

#include "common.h"

#include "BranchPointDB.h"

#include "MBranchPoint.h"
#include "CBranchPoint.h"
#include "Movable.h"
#include "Reactable.h"

//FORWARD DECLARATIONS
class Compartment;
class Cylinder;

/// A container to store a MBranchPoint and CBranchPoint.
/*!
 * BranchPoint class is used to manage and store a MBranchPoint and CBranchPoint.
 * Upon intialization, both of these components are created. Extending the Movable and Reactable
 * classes, the BranchPoint can update its position and reactions according to mechanical equilibration.
 */
class BranchPoint : public Movable, Reactable {
    
private:
    unique_ptr<MBranchPoint> _mBranchPoint; ///< Pointer to mech branch point
    unique_ptr<CBranchPoint> _cBranchPoint; ///< Pointer to chem branch point
    
    Cylinder* _c1; ///< Mother cylinder
    Cylinder* _c2; ///< Branching cylinder
    
    double _position; ///< Position on mother cylinder
    
    short _branchType; ///< Integer specifying the type
    int _branchID; ///< Integer ID of this specific branch point, managed by BranchPointDB
    
    float _birthTime; ///Birth time
    
    Compartment* _compartment; ///< Where this linker is
    
public:
    vector<double> coordinate; ///< coordinate of midpoint, updated with updatePosition()
    
    BranchPoint(Cylinder* c1, Cylinder* c2, short branchType, double position = 0.5, bool creation = false);
    ~BranchPoint();
    
    //@{
    ///Get attached cylinder
    Cylinder* getFirstCylinder() {return _c1;}
    Cylinder* getSecondCylinder() {return _c2;}
    //@}
    
    /// Set chem branch point
    void setCBranchPoint(CBranchPoint* cBranchPoint) {_cBranchPoint = unique_ptr<CBranchPoint>(cBranchPoint);}
    /// Get chem branch point
    CBranchPoint* getCBranchPoint() {return _cBranchPoint.get();}
    
    /// Get mech branch point
    MBranchPoint* getMBranchPoint() {return _mBranchPoint.get();}
    
    //@{
    /// Position management
    double getPosition() {return _position;}
    void setPosition(double position) {_position = position;}
    //@}
    
    //@{
    /// Get br parameter
    short getBranchType() {return _branchType;}
    int getBranchID() {return _branchID;}
    //@}
    
    /// Update the position
    /// @note - changes compartment if needed
    virtual void updatePosition();
    
    /// Update the reaction rates
    virtual void updateReactionRates();
};

#endif
