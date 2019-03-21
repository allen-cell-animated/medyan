
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.2.1
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
#include "CUDAcommon.h"
#include "Bead.h"

//FORWARD DECLARATIONS
class Filament;
class Compartment;
class Bin;

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
                                   public Database< Cylinder > {
    
friend class CController;
friend class DRController;
    
private:
    Bead* _b1;  ///< Pointer to the first bead.
    Bead* _b2; ///< Pointer to the end bead.
    
    unique_ptr<MCylinder> _mCylinder; ///< Pointer to mech cylinder
    unique_ptr<CCylinder> _cCylinder; ///< Pointer to chem cylinder
    
    int _position;          ///< Position on structure
    
    bool _plusEnd = false;  ///< If the cylinder is at the plus end
    bool _minusEnd = false; ///< If the cylinder is at the minus end
    
    short _type; ///< Type of cylinder, either corresponding to Filament or other
                                       
    Compartment* _compartment = nullptr; ///< Where this cylinder is
    
    Cylinder* _branchingCylinder = nullptr; ///< ptr to a branching cylinder
    
    ///For dynamic polymerization rate
    static vector<FilamentRateChanger*> _polyChanger;
                                       
    static ChemManager* _chemManager; ///< A pointer to the ChemManager,
                                      ///< intiailized by CController
    
    ///Helper to get coordinate
    void updateCoordinate();
    
public:
    vector<double> coordinate;
    vector<Bin*> _binvec; //vector of bins. binID corresponding to each binGrid.
    ///< Coordinates of midpoint, updated with updatePosition()
    vector<Bin*> _hbinvec;
    long _dcIndex; ///<Position based on how they occur in Compartment _cylinder vector.
///< Continuous ID assigned for
///< CUDANL calculation
    /// Constructor, initializes a cylinder
    Cylinder(Composite* parent, Bead* b1, Bead* b2, short type, int position,
             bool extensionFront = false,
             bool extensionBack  = false,
             bool initialization = false);
                                       
    virtual ~Cylinder() noexcept;
    
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
    Compartment* getCompartment() {return _compartment;}
    
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
    virtual void addToSubSystem() {}
    virtual void removeFromSubSystem() {
        /* Haoran 03/17/2019
        //Remove from cylinder structure by resetting to default value
        //Reset in bead coordinate vector and add _dbIndex to the list of removedcindex.
        removedcindex.push_back(_dcIndex);
        resetarrays();
        _dcIndex = -1;
        _cylinders.removeElement(this);
        Ncyl = _cylinders.getElements().size();
        */
    }
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
                                       
    virtual void printSelf();
                                       
    /// Returns whether a cylinder is within a certain distance from another
    /// Uses the closest point between the two cylinders
    virtual bool within(Cylinder* other, double dist);

    //Vectorize beads so the coordinates are all available in a single array.
    //@{
    static int maxcindex;
    static int vectormaxsize;
    static int Ncyl; // Currently the value is always numCylinders() - 1
    static vector<int> removedcindex;
    static void revectorizeifneeded(){
        int newsize = vectormaxsize;
        bool check = false;
        if(Bead::triggercylindervectorization || vectormaxsize - maxcindex <= bead_cache/20){

            newsize = (int(Ncyl/cylinder_cache)+2)*cylinder_cache;
            if(removedcindex.size() >= bead_cache)
                newsize = (int(Ncyl/cylinder_cache)+1)*cylinder_cache;
            if(newsize != vectormaxsize || Bead::triggercylindervectorization){
                check = true;
                cylinder* cylindervec = CUDAcommon::serlvars.cylindervec;
                Cylinder** cylinderpointervec = CUDAcommon::serlvars.cylinderpointervec;
                CCylinder** ccylindervec = CUDAcommon::serlvars.ccylindervec;
                delete[] cylindervec;
                delete[] cylinderpointervec;
                delete[] ccylindervec;
                cylinder *newcylindervec = new cylinder[newsize];
                Cylinder **newcylinderpointervec = new Cylinder*[newsize];
                CCylinder **newccylindervec = new CCylinder*[newsize];
                CUDAcommon::serlvars.cylindervec = newcylindervec;
                CUDAcommon::serlvars.cylinderpointervec = newcylinderpointervec;
                CUDAcommon::serlvars.ccylindervec = newccylindervec;
                revectorize(newcylindervec, newcylinderpointervec, newccylindervec);
                vectormaxsize = newsize;
            }
        }
        Bead::triggercylindervectorization = false;
        //@{ check begins
        /*if(check) {
            cylinder *cylindervec = CUDAcommon::serlvars.cylindervec;
            Cylinder **Cylinderpointervec = CUDAcommon::serlvars.cylinderpointervec;
            CCylinder **ccylindervec = CUDAcommon::serlvars.ccylindervec;
            double *coord = CUDAcommon::serlvars.coord;
            std::cout<<"revectorized cylinders"<<endl;
            std::cout << "3 Total Cylinders " << Cylinder::getCylinders().size() << " "
                    "Beads "
                    "" << Bead::getBeads().size() << endl;
            for (auto cyl:Cylinder::getCylinders()) {
                int i = cyl->_dcIndex;
                int id1 = cylindervec[i].ID;
                int id2 = Cylinderpointervec[i]->getID();
                int id3 = ccylindervec[i]->getCylinder()->getID();
                if (id1 != id2 || id2 != id3 || id3 != id1)
                    std::cout << id1 << " " << id2 << " " << id3 << endl;
                auto b1 = cyl->getFirstBead();
                auto b2 = cyl->getSecondBead();
                long idx1 = b1->getIndex();
                long idx2 = b2->getIndex();
                cylinder c = cylindervec[i];
                std::cout << "3 bindices for cyl with ID "<<cyl->getID()<<" cindex " << i <<
                " are "<< idx1 << " " << idx2 << " " << c.bindices[0] << " " << c.bindices[1] << endl;
                if (c.bindices[0] != idx1 || c.bindices[1] != idx2) {

                    std::cout << "Bead " << b1->coordinate[0] << " " << b1->coordinate[1]
                              << " " << b1->coordinate[2] << " " << " " << b2->coordinate[0]
                              << " " << b2->coordinate[1] << " " << b2->coordinate[2]
                              << " idx " << b1->getIndex() << " " << b2->getIndex() << endl;

                    std::cout << coord[3 * idx1] << " " << coord[3 * idx1 + 1] << " "
                              << coord[3 * idx1 + 2] << " " << coord[3 * idx2] << " "
                              << coord[3 * idx2 + 1] << " " << coord[3 * idx2 + 2] << endl;
                }
            }
        }*/
        //@} check ends.
    }
    static void revectorize(cylinder* cylindervec, Cylinder** cylinderpointervec,
                            CCylinder** ccylindervec);
    void  copytoarrays();
    void resetarrays();
    void resetcylinderstruct(cylinder* cylindervec, long idx){
        cylindervec[idx].filamentID = -1;
        cylindervec[idx].filamentposition = -1;
        cylindervec[idx].beads[0]= nullptr;
        cylindervec[idx].beads[1]= nullptr;
        cylindervec[idx].cmpID = -1;
        cylindervec[idx].cindex = -1;
        cylindervec[idx].coord[0] = -1.0;
        cylindervec[idx].coord[1] = -1.0;
        cylindervec[idx].coord[2] = -1.0;
        cylindervec[idx].type = -1;
        cylindervec[idx].ID = -1;
        cylindervec[idx].availbscount = -1;
    }
    //@}
};

#endif
