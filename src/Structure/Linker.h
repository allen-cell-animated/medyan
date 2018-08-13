
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_Linker_h
#define MEDYAN_Linker_h

#include "common.h"

#include "Composite.h"
#include "CLinker.h"
#include "MLinker.h"

#include "Database.h"
#include "Histogram.h"
#include "Trackable.h"
#include "Movable.h"
#include "Reactable.h"
#include "RateChangerImpl.h"

//FORWARD DECLARATIONS
class Cylinder;
class Controller;
class DRController;

/// A container to store a MLinker and CLinker.
/*!
 *  Linker class is used to manage and store a MLinker and CLinker. Upon intialization,
 *  both of these components are created.
 *
 *  Extending the Movable class, the positions of all instances 
 *  can be updated by the SubSystem.
 *
 *  Extending the Reactable class, the reactions associated with all
 *  instances can be updated by the SubSystem.
 */
class Linker : public Component, public Trackable, public Movable, public Reactable {

friend class Controller;
friend class DRController;
    
private:
    unique_ptr<MLinker> _mLinker; ///< Pointer to mech linker
    unique_ptr<CLinker> _cLinker; ///< Pointer to chem linker
    
    Cylinder* _c1; ///< First cylinder the linker is bound to
    Cylinder* _c2; ///< Second cylinder the linker is bound to
    
    double _position1; ///< Position on first cylinder
    double _position2; ///< Position on second cylinder
    
    short _linkerType; ///< Integer specifying the type
    int _linkerID; ///< Integer ID of this specific linker, managed by Database
    
    float _birthTime; ///Birth time
    
    Compartment* _compartment; ///< Where this linker is
    
    static Database<Linker*> _linkers;
    ///< Collection in SubSystem
    
    //@{
    ///Histogram data
    static Histogram* _lifetimes;
    //@}
    
    ///For dynamic rate unbinding
    static vector<LinkerRateChanger*> _unbindingChangers;
    
    ///Helper to get coordinate
    void updateCoordinate();
    
public:
    vector<double> coordinate;
    ///< coordinate of midpoint, updated with updatePosition()

    
    Linker(Cylinder* c1, Cylinder* c2, short linkerType,
           double position1 = 0.5, double position2 = 0.5);
    
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
    virtual int getType() {return _linkerType;}
    int getID() {return _linkerID;}
    //@}
    
    /// Get the birth time
    float getBirthTime() {return _birthTime;}
    
    //@{
    /// SubSystem management, inherited from Trackable
    virtual void addToSubSystem() { _linkers.addElement(this);}
    virtual void removeFromSubSystem() {_linkers.removeElement(this);}
    //@}
    
    /// Get all instances of this class from the SubSystem
    static const vector<Linker*>& getLinkers() {
        return _linkers.getElements();
    }
    /// Get the number of linkers in this system
    static int numLinkers() {
        return _linkers.countElements();
    }
    
    /// Get the lifetimes
    static Histogram* getLifetimes() {return _lifetimes;}
    
    /// Update the position, inherited from Movable
    /// @note - changes compartment if needed
    virtual void updatePosition();
    
    /// Update the reaction rates, inherited from Reactable
    virtual void updateReactionRates();
    
    virtual void printSelf();
    
    /// Count the number of linker species with a given name in the system
    static species_copy_t countSpecies(const string& name);
};


#endif 
