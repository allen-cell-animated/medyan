
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

#ifndef M3SYM_Utility_h
#define M3SYM_Utility_h

#include <memory>

using namespace std;

template<typename T, typename ...Args>
unique_ptr<T> make_unique( Args&& ...args )
{
    return unique_ptr<T>( new T( forward<Args>(args)... ) );
}

template<typename T, typename U>
bool isSame(const U& x) {
    return typeid(x) == typeid(T&); 
}

#endif
