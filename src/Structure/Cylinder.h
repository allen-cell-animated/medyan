
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_Cylinder_h
#define MEDYAN_Cylinder_h

#include <iostream>

#include "common.h"

#include "MCylinder.h"
#include "CCylinder.h"
#include "RateChanger.h"

#include "Database.h"
#include "Trackable.h"
#include "Movable.h"
#include "Reactable.h"
#include "DynamicNeighbor.h"
#include "Component.h"
#include "Bead.h"
#include "Structure/CellList.hpp"
#include "Util/Math/Vec.hpp"

//FORWARD DECLARATIONS
class Filament;
class Compartment;
class Bin;

struct CylinderInfoData {
    struct CylinderInfo {
        int filamentId = -1;
        int positionOnFilament = -1;
        int filamentFirstEntry = 0;
        int compartmentId = -1;
        std::size_t beadIndices[2];
        mathfunc::Vec< 3, floatingpoint > coord;
        short type = -1;
        int id = -1;
        CCylinder* chemCylinder;
    };

    std::vector< CylinderInfo > value;

    void push_back(CylinderInfo c) { value.push_back(c); }
    void set_content(std::size_t pos, CylinderInfo c) { value[pos] = c; }
    void move_content(std::size_t from, std::size_t to) { value[to] = value[from]; }
    void resize(std::size_t size) { value.resize(size); }
};

/// A container to store a MCylinder and CCylinder.
/*!
 *  Cylinder class is used to manage and store a MCylinder and CCylinder.
 *  Upon intialization, both of these components are created.
 *
 *  Extending the Movable class, the positions of all instances
 *  can be updated by the SubSystem.
 *
 *  Extending the Reactable class, the reactions associated with
 *  all instances can be updated by the SubSystem.
 *
 *  Extending the DynamicNeighbor class, all instances can be
 *  kept in [NeighborLists](@ref NeighborList).
 */
class Cylinder : public Component, public Trackable, public Movable,
                                   public Reactable, public DynamicNeighbor,
                                   public Database< Cylinder, true, CylinderInfoData > {

friend class CController;
friend class DRController;

private:

    chrono::high_resolution_clock::time_point mins, mine;

    Bead* _b1;  ///< Pointer to the first bead.
    Bead* _b2; ///< Pointer to the end bead.

    unique_ptr<MCylinder> _mCylinder; ///< Pointer to mech cylinder
    unique_ptr<CCylinder> _cCylinder; ///< Pointer to chem cylinder

    bool _plusEnd = false;  ///< If the cylinder is at the plus end
    bool _minusEnd = false; ///< If the cylinder is at the minus end

    short _type; ///< Type of cylinder, either corresponding to Filament or other

	int _position;          ///< Position on structure

    cell_list::CellListElementUser< Cylinder, Compartment > _cellElement;

    Cylinder* _branchingCylinder = nullptr; ///< ptr to a branching cylinder

    ///For dynamic polymerization rate
    static vector<FilamentRateChanger*> _polyChanger;

    static ChemManager* _chemManager; ///< A pointer to the ChemManager,
                                      ///< intiailized by CController

    ///Helper to get coordinate
    void updateCoordinate();


    /// ID of filament
    int _filID;

public:
    using DatabaseType = Database< Cylinder, true, CylinderInfoData >;

    vector<floatingpoint> coordinate;
    vector<Bin*> _binvec; //vector of bins. binID corresponding to each binGrid.
    ///< Coordinates of midpoint, updated with updatePosition()
    vector<Bin*> _hbinvec;

    static bool setpositionupdatedstate; //Setter to check if position has been updated

    // Update CylinderInfoData using newest information in the system
    static void updateAllData() {
        // Update data for all cylinders
        for(auto c : getCylinders()) c->updateData();
    }
    void updateData(); // Update data for this cylinder. TODO: make it const

    /// Constructor, initializes a cylinder
    Cylinder(Composite* parent, Bead* b1, Bead* b2, short type, int position,
             bool extensionFront = false,
             bool extensionBack  = false,
             bool initialization = false);

    virtual ~Cylinder() noexcept;

    const auto& getCoordinate() const { return getDbData().value[getStableIndex()].coord; }
    auto      & getCoordinate()       { return getDbData().value[getStableIndex()].coord; }

    /// Get mech cylinder
    MCylinder* getMCylinder() {return _mCylinder.get();}

    /// Get chem cylinder
    CCylinder* getCCylinder() {return _cCylinder.get();}
    /// set chem cylinder
    /// @note: since this is a unique ptr, will implicitly delete old chem cylinder
    void setCCylinder(CCylinder* c) {_cCylinder = unique_ptr<CCylinder>(c);}

    /// Get cylinder type
    virtual int getType();

    //@{
    /// Get beads
    Bead* getFirstBead() {return _b1;}
    Bead* getSecondBead() {return _b2;}
    //@}

    //@{
    /// Set beads
    void setFirstBead(Bead* b) {_b1 = b;}
    void setSecondBead(Bead* b) {_b2 = b;}
    //@}

    /// Get compartment
    Compartment* getCompartment() const { return _cellElement.manager->getHeadPtr(_cellElement); }

    //@{
    /// Branching cylinder management
    Cylinder* getBranchingCylinder() {return _branchingCylinder;}
    void setBranchingCylinder(Cylinder* c) {_branchingCylinder = c;}
    //@}

    ///@{
    /// Set plus and minus end boolean markers
    bool isPlusEnd() {return _plusEnd;}
    void setPlusEnd(bool plusEnd) {_plusEnd = plusEnd;}

    bool isMinusEnd() {return _minusEnd;}
    void setMinusEnd(bool minusEnd) {_minusEnd = minusEnd;}
    //@}

    int getPosition() {return _position;}

    //@{
    /// SubSystem management, inherited from Trackable
    // Does nothing
    virtual void addToSubSystem() override {}
    virtual void removeFromSubSystem() override {}
    //@}

    /// Get all instances of this class from the SubSystem
    static const vector<Cylinder*>& getCylinders() {
        return getElements();
    }
    /// Get the number of cylinders in this system
    static int numCylinders() {
        return getElements().size();
    }

    /// Update the position, inherited from Movable
    /// @note - changes compartment if needed
    virtual void updatePosition();

    /// Update the reaction rates, inherited from Reactable
    virtual void updateReactionRates();

    /// Check if this cylinder is grown to full length
    bool isFullLength();

    virtual void printSelf()const;

    /// Returns whether a cylinder is within a certain distance from another
    /// Uses the closest point between the two cylinders
    virtual bool within(Cylinder* other, floatingpoint dist);
    void setFilID(int filID){
        _filID = filID;
    };

    int getFilID() const {
        return _filID;
    }

    static floatingpoint timecylinder1;
	static floatingpoint timecylinder2;
	static floatingpoint timecylinderchem;
	static floatingpoint timecylindermech;
    //@}
};

#endif
