
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

#ifndef MEDYAN_utility_h
#define MEDYAN_utility_h

#include <tuple>
#include <chrono>
#include <memory>
#include <math.h>
#include <vector>
#include <sstream>
#include <iterator>

#ifdef CUDAACCL
#include <cuda.h>
#include <cuda_runtime.h>
#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#else
#define CUDA_HOSTDEV
#endif

#endif
using namespace std;

///floatingpoint typedef
typedef float floatingpoint;
typedef double doubleprecision;

//to test for zero values
const floatingpoint ZERO_PREC = 1E-4;

/// A random seed based on clock cycles
extern unsigned long long rdtsc();

///Check equaility of floatingpoints
#ifdef CUDAACCL
__host__ __device__
#endif
inline bool areEqual(floatingpoint d1, floatingpoint d2) {
    const floatingpoint ZERO_PREC = 1E-4;
    return fabs(d1 - d2) < ZERO_PREC;
}

//checking for equality with a lower threshold.
inline bool areEqualLT(floatingpoint d1, floatingpoint d2) {
	const floatingpoint ZERO_PRECLT = 1E-1;
	return fabs(d1 - d2) < ZERO_PRECLT;
}

/// Safe arc cos function
#ifdef CUDAACCL
 __host__ __device__
#endif
inline floatingpoint safeacos (floatingpoint x) {
    if (x < -1.0) x = -1.0;
    else if (x > 1.0) x = 1.0;
    return acos(x);
}

/// Make a unique ptr
#if __cplusplus <= 201103L
template<typename T, typename ...Args>
unique_ptr<T> make_unique( Args&& ...args )
{
    return unique_ptr<T>( new T( forward<Args>(args)... ) );
}
#endif
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
