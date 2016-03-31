
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2015)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef M3SYM_utility_h
#define M3SYM_utility_h

#include <tuple>
#include <chrono>
#include <memory>
#include <math.h>
#include <vector>
#include <sstream>
#include <iterator>

using namespace std;

//to test for zero values
const double ZERO_PREC = 1E-6;

/// A random seed based on clock cycles
extern unsigned long long rdtsc();

///Check equaility of doubles
inline bool areEqual(double d1, double d2) {
    
    return fabs(d1 - d2) < ZERO_PREC;
}

/// Safe arc cos function
inline double safeacos (double x) {
    if (x < -1.0) x = -1.0;
    else if (x > 1.0) x = 1.0;
    return acos(x);
}

/// Make a unique ptr
template<typename T, typename ...Args>
unique_ptr<T> make_unique( Args&& ...args )
{
    return unique_ptr<T>( new T( forward<Args>(args)... ) );
}

/// Compare types
template<typename T, typename U>
bool isSame(const U& x) {
    return typeid(x) == typeid(T&); 
}

/// Split a string by whitespace into generic type
template<typename T>
vector<T> split(const string& line) {
    istringstream is(line);
    return vector<T>(istream_iterator<T>(is), istream_iterator<T>());
}

//Functions to hash by a tuple
namespace std{
    namespace
    {
        // Code from boost
        // Reciprocal of the golden ratio helps spread entropy
        // and handles duplicates.
        // See Mike Seymour in magic-numbers-in-boosthash-combine:
        // http://stackoverflow.com/questions/4948780
        
        template <class T>
        inline void hash_combine(std::size_t& seed, T const& v)
        {
            seed ^= hash<T>()(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }
        
        // Recursive template code derived from Matthieu M.
        template <class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
        struct HashValueImpl
        {
            static void apply(size_t& seed, Tuple const& tuple)
            {
                HashValueImpl<Tuple, Index-1>::apply(seed, tuple);
                hash_combine(seed, get<Index>(tuple));
            }
        };
        
        template <class Tuple>
        struct HashValueImpl<Tuple,0>
        {
            static void apply(size_t& seed, Tuple const& tuple)
            {
                hash_combine(seed, get<0>(tuple));
            }
        };
    }
    
    template <typename ... TT>
    struct hash<std::tuple<TT...>>
    {
        size_t
        operator()(std::tuple<TT...> const& tt) const
        {
            size_t seed = 0;
            HashValueImpl<std::tuple<TT...> >::apply(seed, tt);
            return seed;
        }
        
    };
    
    ///Sum a vector of shorts
    inline short sum(vector<short> vec) {
        
        short sum = 0;
        for (auto val : vec) sum += val;
        
        return sum;
    }
}

#endif
