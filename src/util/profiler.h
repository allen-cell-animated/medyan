#ifndef MEDYAN_UTIL_PROFILER_H
#define MEDYAN_UTIL_PROFILER_H

#include <chrono>
#include <ostream>
#include <string>
#include <type_traits>

namespace profiler {

#ifdef USE_PROFILER
constexpr bool useProfiler = true;
#else
constexpr bool useProfiler = false;
#endif

namespace internal {

template< bool > class SimpleTimerImpl;

// Specialization
template<> class SimpleTimerImpl< true > {
private:
    using _time_point_t = std::chrono::time_point< std::chrono::steady_clock >;
    using _time_diff_t = decltype(_time_point_t() - _time_point_t());

    _time_point_t _start;
    _time_diff_t _elapsed;

    std::string _name;
public:
    // Record the current time as start time.
    void start() {
        _start = std::chrono::steady_clock::now();
    }
    // Calculate the time elapsed from start to current time.
    void elapse() {
        _elapsed = std::chrono::steady_clock::now() - _start;
    }
    void report();
};

template<> class SimpleTimerImpl< false > {
public:
    void start() {}
    void elapse() {}
    void report() {}
};

} // namespace internal

using Timer = internal::SimpleTimerImpl< useProfiler >;

} // namespace profiler

#endif
