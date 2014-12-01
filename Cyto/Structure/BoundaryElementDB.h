
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

#ifndef M3SYM_BoundaryElementDB_h
#define M3SYM_BoundaryElementDB_h

#include "common.h"

#include "BoundaryElementImpl.h"

/// BoundaryElementDB class is a database for all [BoundaryElements](@ref BoundaryElement) in the system
/*!
 *   This BoundaryElementDB inherits from list and manage all creations and removing of
 *   [BoundaryElement](@ref BoundaryElement) objects, as well as some standard list functions and iterators.
 *   The [BoundarySurface] (@ref BoundarySurface) class calls this database to create and/or remove boundary elements.
 */
class BoundaryElementDB: private list<BoundaryElement*>
{
    typedef list<BoundaryElement*> bedb;
    
public:
    using bedb::size;
    using bedb::begin;
    using bedb::end;
    using bedb::erase;
    using bedb::remove;
    
    /// Copying is not allowed
    BoundaryElementDB(const BoundaryElementDB &rhs) = delete;
    
    /// Assignment is not allowed
    BoundaryElementDB& operator=(BoundaryElementDB &rhs) = delete;
    
    /// Get the instance of this singleton
    static BoundaryElementDB* instance();
    
    /// Create a new plane boundary element
    BoundaryElement* createPlaneBoundaryElement(vector<double>& coords, vector<double>& normal,
                                                double repulsConst, double screenLength) {
        
        BoundaryElement* b = new PlaneBoundaryElement(coords, normal, repulsConst, screenLength);
        push_back(b);
        return b ;
    }
    
    /// Create a spherical boundary element
    BoundaryElement* createSphereBoundaryElement(vector<double>& coords, double radius,
                                                 double repulsConst, double screenLength) {
        
        BoundaryElement* b = new SphereBoundaryElement(coords, radius, repulsConst, screenLength);
        push_back(b);
        return b ;
    }
    
    /// Create a cylindrical z boundary element
    BoundaryElement* createCylindricalZBoundaryElement(vector<double> coords, double radius,
                                                       double height, double repulsConst, double screenLength) {
        
        BoundaryElement* b = new CylindricalZBoundaryElement(coords, radius, height, repulsConst, screenLength);
        push_back(b);
        return b ;
    }
    
    /// Create a half sphere z boundary element
    BoundaryElement* createHalfSphereZBoundaryElement(vector<double> coords, double radius,
                                                      bool up, double repulsConst, double screenLength) {
        
        BoundaryElement* b = new HalfSphereZBoundaryElement(coords, radius, up, repulsConst, screenLength);
        push_back(b);
        return b ;
    }
    
    /// Remove boundary element
    void removeBoundaryElement(BoundaryElement* b){
        
        remove(b);
        delete b;
        
    }
private:
    static BoundaryElementDB* _instance; ///< Singleton instance
    BoundaryElementDB() {};
    
};








#endif /* defined(__Cyto__BoundaryElementDB__) */
