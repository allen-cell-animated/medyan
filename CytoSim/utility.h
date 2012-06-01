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

template<typename T1>
struct isSame {
    template<typename T2>
    bool operator() (const T2& t2, typename std::enable_if<std::is_same<T1, T2>::value>::type* p = nullptr){
        return true;
    }
    
    template<typename T2>
    bool operator() (const T2& t2, typename std::enable_if<!std::is_same<T1, T2>::value>::type* p = nullptr) 
    { 
        return false;
    }                                                            
};

#endif
