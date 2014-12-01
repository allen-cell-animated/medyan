
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

#ifndef M3SYM_FilamentDB_h
#define M3SYM_FilamentDB_h

#include <list>

#include "common.h"

#include "Filament.h"

#include "MathFunctions.h"
#include "SystemParameters.h"


/// FilamentDB class is a database for all boundary elements in the system
/*!
 *   This FilamentDB inherits from list and manage all creations and removing of
 *   [Filaments](@ref Filament) objects, as well as some standard list functions and iterators.
 *   The [SubSystem] (@ref SubSystem) class calls this database to create and/or remove filaments.
 */
class FilamentDB: private list<Filament*> {
    typedef list<Filament*> fdb;
    
public:
    using fdb::size;
    using fdb::begin;
    using fdb::end;
    
    static FilamentDB* instance();
    
    /// Create a filament, given a vector of initial bead coordinates
    Filament* createFilament(SubSystem* s, vector<vector<double> >& v) {
        
        double d = mathfunc::TwoPointDistance(v[0], v[1]);
        vector<double> tau = mathfunc::TwoPointDirection(v[0], v[1]);
        
        int numSegment = d / SystemParameters::Geometry().cylinderSize;
        
        // check how many segments can fit between end-to-end of the filament
        if (numSegment == 0){
            
            Filament* pf = new Filament(s, v[0], tau, _currentFilamentID++); //create a filament with only two beads
            push_back(pf);

            return pf;
        }
        
        else {
            Filament* pf = new Filament(s, v, numSegment + 1, _currentFilamentID++, "STRAIGHT");  //Create a long filament with numSeg.
            push_back(pf);
            
            return pf;
        }
    }

    /// Remove a filament from the system
    void removeFilament(Filament* f) {
        delete f;
        remove(f);
    };
    
private:
    static int _currentFilamentID;  ///< To assign filament ids

    static FilamentDB* _instance;   ///< Singleton instance
    FilamentDB() {};
    
};

#endif
