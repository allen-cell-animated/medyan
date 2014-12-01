
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

#ifndef M3SYM_BeadDB_h
#define M3SYM_BeadDB_h

#include <list>
#include <vector>

#include "common.h"

#include "Bead.h"

/// BeadDB class is a database for all [Beads](@ref Bead) in the system
/*!  
 *   This BeadDB inherits from list and manage all creations and removing of
 *   [Bead](@ref Bead) objects, as well as some standard list functions and iterators.
 *   The [Filament] (@ref Filament) class calls this database to create and/or remove beads.
 */
class BeadDB: private list<Bead*>
{
    typedef list<Bead*> bdb;
    
public:
    using bdb::size;
    using bdb::begin;
    using bdb::end;
    using bdb::erase;
    using bdb::remove;
    
    /// Copying is not allowed
    BeadDB(const BeadDB &rhs) = delete;
    
    /// Assignment is not allowed
    BeadDB& operator=(BeadDB &rhs) = delete;
    
    /// Get the instance of this singleton
    static BeadDB* instance();
    
    /// Create a new bead with no coordinates
    Bead* createBead(int positionFilament) {
        
        Bead* b = new Bead(positionFilament);
        push_back(b);
        return b ;}
    
    /// Create bead with a given coordinate on a given filament
    Bead* createBead(vector<double>& v, int positionFilament) {
        
        Bead* b = new Bead(v, positionFilament);
        push_back(b);
        return b ;}
    
    
    /// Remove bead
    void removeBead(Bead* b){
        delete b;
        remove(b);
        
    }
private:
    static BeadDB* _instance; ///< Singleton instance
    BeadDB() {};
};

#endif 
