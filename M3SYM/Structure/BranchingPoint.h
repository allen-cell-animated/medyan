
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

#ifndef M3SYM_BranchingPoint_h
#define M3SYM_BranchingPoint_h

#include "common.h"

#include "BranchingPointDB.h"

#include "MBranchingPoint.h"
#include "CBranchingPoint.h"
#include "Movable.h"

//FORWARD DECLARATIONS
class Compartment;
class Cylinder;

/// A container to store a MBranchingPoint and CBranchingPoint.
/*!
 * BranchingPoint class is used to manage and store a MBranchingPoint and 
 * CBranchingPoint. Upon intialization, both of these components are created. Extending
 * the Movable class, the BranchingPoint can update its position according to mechanical 
 * equilibration.
 */
class BranchingPoint : public Movable{
    
private:
    unique_ptr<MBranchingPoint> _mBranchingPoint; ///< Pointer to mech branch point
    unique_ptr<CBranchingPoint> _cBranchingPoint; ///< Pointer to chem branch point
    
    Cylinder* _c1; ///< Mother cylinder
    Cylinder* _c2; ///< Branching cylinder
    
    double _position; ///< Position on mother cylinder
    
    short _branchType; ///< Integer specifying the type
    int _branchID; ///< Integer ID of this specific branch point,
                   ///< managed by BranchingPointDB
    
    float _birthTime; ///Birth time
    
    Compartment* _compartment; ///< Where this linker is
    
public:
    vector<double> coordinate; ///< coordinate of midpoint,
                               ///< updated with updatePosition()
    
    BranchingPoint(Cylinder* c1, Cylinder* c2,
                   short branchType, double position = 0.5,
                   bool creation = false);
    ~BranchingPoint();
    
    //@{
    ///Get attached cylinder
    Cylinder* getFirstCylinder() {return _c1;}
    Cylinder* getSecondCylinder() {return _c2;}
    //@}
    
    /// Set chem branch point
    void setCBranchingPoint(CBranchingPoint* cBranchingPoint) {
        _cBranchingPoint = unique_ptr<CBranchingPoint>(cBranchingPoint);
    }
    /// Get chem branch point
    CBranchingPoint* getCBranchingPoint() {return _cBranchingPoint.get();}
    
    /// Get mech branch point
    MBranchingPoint* getMBranchingPoint() {return _mBranchingPoint.get();}
    
    //@{
    /// Position management
    double getPosition() {return _position;}
    void setPosition(double position) {_position = position;}
    //@}
    
    //@{
    /// Get branch parameter
    short getBranchType() {return _branchType;}
    int getBranchID() {return _branchID;}
    //@}
    
    /// Update the position
    /// @note - changes compartment if needed
    virtual void updatePosition();
};

#endif
