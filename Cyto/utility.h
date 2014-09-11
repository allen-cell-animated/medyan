//
//  utility.h
//  CytoSim
//
//  Created by Garegin Papoian on 5/4/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_utility_h
#define CytoSim_utility_h

#include <memory>

template<typename T, typename ...Args>
std::unique_ptr<T> make_unique( Args&& ...args )
{
    return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
}

template<typename T, typename U>
bool isSame(const U& x) {
    return typeid(x) == typeid(T&); 
}

#endif
