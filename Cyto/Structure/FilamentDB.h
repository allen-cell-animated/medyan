//
//  FilamentDB.h
//  CytoMech
//
//  Created by Konstantin Popov on 4/15/14.
//  Copyright (c) 2014 Konstantin Popov. All rights reserved.
//

#ifndef CytoMech_FilamentDB_h
#define CytoMech_FilamentDB_h


#include <iostream>
#include <list>

#include "common.h"
#include "MathFunctions.h"
#include "SystemParameters.h"
#include "Filament.h"

///Key to access instance of FilamentDB
class FilamentDBKey {friend class SubSystem;
                    friend class FilamentFF;
                    friend class MController;
                    friend class Output;
#ifdef TESTING
                    public:
#endif //TESTING
                    FilamentDBKey(){};
                    ~FilamentDBKey(){}; };


///FilamentDB is used to store all filaments in the system
/*! An Object Data Base singleton structure will be used as a container for all main objects: Beads, Filament, Linkers,
 *  Boundary Elements, and Motors. This structure inherits from std:: list and manage all creations and removing 
 *  of objects, as well as some stabdart list functions and iterators.
 */
class FilamentDB: private std::list<Filament*> {
    typedef std::list<Filament*> fdb;
    
public:
    using fdb::size;
    using fdb::begin;
    using fdb::end;
    
    static FilamentDB* Instance(FilamentDBKey k);
    
    Filament* CreateFilament(SubSystem* s, std::vector<std::vector<double> >& v) {
        
        double d = mathfunc::TwoPointDistance(v[0], v[1]);
        std::vector<double> tau = mathfunc::TwoPointDirection(v[0], v[1]);
        
        int numSegment = d/SystemParameters::Geometry().cylinderSize;
        // check how many segments can fit between end-to-end of the filament
        
        if (numSegment == 0){
            
            Filament* pf = new Filament(s, v[0], tau, _currentFilamentID++); //create a filament with only two beads
            push_back(pf);
            std::cout<<"short filament created"<<std::endl;
            return pf;}
        
        else {
            Filament* pf = new Filament(s, v, numSegment + 1, _currentFilamentID++, "STRAIGHT");  //Create a long filament with numSeg.
            push_back(pf);
            std::cout<<"long filament created"<<std::endl;
            
            return pf;
        }
    }

    void RemoveFilament(Filament* f) {
        delete f;
        remove(f);
    };
    
private:
    ///To assign filament ids
    static int _currentFilamentID;

    static FilamentDB* _instance;
    FilamentDB() {};
    
};

#endif
