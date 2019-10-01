#ifndef MEDYAN_Util_ThreadUtil_hpp
#define MEDYAN_Util_ThreadUtil_hpp

#include <algorithm> // min

#include "ThreadPool.hpp"

// This file implements some helper functions that might help with
// thread-based parallelization.
//
// The implementation is based on threadings using the thread pool. The auto
// vectorization is not used in this file.
//
// Some of the ideas (Parallel TS) are present starting C++17. It is possible
// that the thread pooling and this file will no longer be needed with future
// versions of C++ and compilers.

template< typename I, typename F >
inline void poolForkJoinFixedSize(ThreadPool& tp, I begin, I end, I batchSize, F&& f) {
    // F must take the indexing int as the parameter

    using ReturnType = std::result_of_t< F(I) >;

    std::vector< std::future< ReturnType > > futs;
    for(I i = begin; i < end; i += batchSize) {
        I thisBegin = i;
        I thisEnd   = std::min<I>(end, i + batchSize);

        futs.push_back(tp.submit(
            [ f { std::forward<F>(f) } ] (I begin, I end) {
                for(I i = begin; i < end; ++i) f(i);
            },
            thisBegin, thisEnd
        ));
    }

    // wait for finish
    for(auto& f : futs) f.wait();
}

template< typename I, typename F >
inline auto poolForkJoinAccumulateFixedSize(ThreadPool& tp, I begin, I end, I batchSize, F&& f) {
    // F must take the indexing int as the parameter.
    // F must return a type that is associative under "+" operation.

    using ReturnType = std::result_of_t< F(I) >;

    std::vector< std::future< ReturnType > > futs;
    for(I i = begin; i < end; i += batchSize) {
        I thisBegin = i;
        I thisEnd   = std::min<I>(end, i + batchSize);

        futs.push_back(tp.submit(
            [ f { std::forward<F>(f) } ] (I begin, I end) {
                ReturnType res = 0;
                for(I i = begin; i < end; ++i) res += f(i);
                return res;
            },
            thisBegin, thisEnd
        ));
    }

    // accumulate the results
    ReturnType res = 0;
    for(auto& f : futs) res += f.get();
    return res;
}

#endif
