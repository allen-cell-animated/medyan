#ifndef MEDYAN_UTIL_PROFILER_H
#define MEDYAN_UTIL_PROFILER_H

#include <chrono>
#include <string>

namespace profiler {

#ifdef USE_PROFILER
constexpr bool useProfiler = true;
#else
constexpr bool useProfiler = false;
#endif

namespace internal {

template<
    bool enable,
    typename TimerConf
> class TimerImpl {
private:
    std::chrono::time_point< std::chrono::steady_clock > _start;
    std::chrono::time_point< std::chrono::steady_clock > _end;
public:
    // Record the current time as start time.
    void start() {
        _start = std::chrono::steady_clock::now();
    }
};

// Specialization
template< typename TimerConf > class TimerImpl< false, TimerConf > {};

} // namespace internal

using Timer = internal::TimerImpl< useProfiler, void >;

} // namespace profiler

#endif
